#' @title Normalize RNA-seq count table
#' @author Jacques van Helden and Mustafa AbuElqumsan
#' @description Normalize a table of raw counts from RNA-seq
#' @param counts a raw count table, with one column per sample and one row per gene.
#' Normalization will be done by samples (columns).
#' @param phenoTable the pheno table associated to the count table. Should be provided in order to adapt it in case samples would be suppressed because they have a scaling factor of 0.
#' @param classLabels the class labels associated each sample of the count table. Should be provided in order to adapt it in case samples would be suppressed because they have a scaling factor of 0.
#'
#' @param method="quantile" normalization method. Supported: sum, mean, media, quantile.
#' Preferred method: quantile = 0.75. Quantiles are more robust to outliers than the mean
#' or sum. The median is sometimes weak in RNA-seq data, the 75th percentile (quantile 0.75)
#' seems a good tradeoff between robustness and representativity (taking into account a
#' representative proportion of the total counts per samples).
#' Note that median is the same as quantile 0.5.
#'
#' @param quantile=0.75 quantile used as scaling factor when quantile method is selected.
#' @param log2=TRUE apply a log2 transformation
#' @param epsilon=0.1 value added to all counts before applying the log2 transformation
#' in order to avoid zero counts.
#'
#' @return a list with the normalized count table + a table providing sample-wise statistics
#'
#' @examples
#'
#' @export
NormalizeCounts <- function(counts,
                            phenoTable,
                            classLabels,
                            method="quantile",
                            quantile=0.75,
                            log2=TRUE,
                            epsilon=0.1) {
  message.with.time("Starting NormalizeCounts() for Recount experiment ID ", parameters$recountID)

  ## Compute sample-wise statistics
  sampleStats <- data.frame(
    sum = apply(counts, 2, sum, na.rm=TRUE),
    min = apply(counts, 2, min, na.rm=TRUE),
    mean = apply(counts, 2, mean, na.rm=TRUE),
    Q1 = apply(counts, 2, quantile, na.rm=TRUE, probs=0.25),
    median = apply(counts, 2, median, na.rm=TRUE),
    Q3 = apply(counts, 2, quantile, na.rm=TRUE, probs=0.75),
    max = apply(counts, 2, max, na.rm=TRUE),
    zero.values = apply(counts == 0, 2, sum),
    na.values = apply(is.na(counts), 2, sum),
    infinite.values = apply(is.infinite(counts), 2, sum)
  )

  ## check it the original data contains NA values
  if (sum(sampleStats$na.values) > 0) {
    message ("\tRaw counts contain ", sum(sampleStats$na.values), " NA values. ")
  }

  if (method == "quantile") {
    sampleStats$norm.quantile <- apply(counts, 2, quantile, na.rm=TRUE, probs=quantile)
    sampleStats$scaling.factor <- sampleStats$norm.quantile
  } else if (method == "mean") {
    ## Note: mean is very sensitive to outliers, which are very problematic with RNAseq data
  #  sampleStats$scaling.factor <- apply(counts, 2, quantile, na.rm=TRUE, probs=quantile)
    sampleStats$scaling.factor <- sampleStats$mean
  } else if (method == "median") {
    sampleStats$scaling.factor <- sampleStats$median
  } else if (method == "sum") {
    sampleStats$scaling.factor <- sampleStats$sum
  } else {
    stop(method, " is not a valid method for NormalizeCounts()")
  }

  ## Detect problems related to null scaling factors, which may happen in some datasets due to a very large number of zeros.
  zeroScaledSamples <- sampleStats$scaling.factor == 0
  discaredSampleNames <- colnames(counts[, zeroScaledSamples])

  # apply(counts[, zeroScaledSamples]==0, 2, sum)
  message("\tDiscarding ", sum(zeroScaledSamples), " samples because their scaling factor is null. ")
  message("\tDiscarded samples: ", paste(collapse="; ", discaredSampleNames))

  ## Compute normalised counts
  normTarget <- median(sampleStats$scaling.factor[!zeroScaledSamples])
  normCounts <- t(t(counts[, !zeroScaledSamples]) / sampleStats$scaling.factor[!zeroScaledSamples]) * normTarget
  # Check that all normalised coutns have the same scaling factor
  # apply(normCounts, 2, quantile, probs=quantile)

  ## Run log2 transformation if required
  if (log2) {
    ## TO BE CHANGED
    ## In principle there is no reason for the log2 transformation to generate NA values.
    ## We should certainly not use na.omit, because this removes any sample that would have any NA value.
    ## Instead, we need to understand why there are NA values in the original count table.
    normCounts <- log2(normCounts + epsilon)
    # dim(normCounts)
    # normCounts <-na.omit(normCounts) ## THIS WAS REDOING THE SAME AS ABOVE

  }


  # hist(unlist(normCounts), breaks=100)

  ## Build the result list
  result <- list()
  result$method <- method
  if (method == "quantile") {
    result$quantile <- quantile
  }
  result$sampleStats <- sampleStats
  result$discaredSampleNames <- discaredSampleNames

  result$nb.samples <- ncol(normCounts)
  result$nb.genes <- nrow(normCounts)
  result$counts <- normCounts
  result$phenoTable <- phenoTable[!zeroScaledSamples,]
  result$classLabels <- classLabels[!zeroScaledSamples]


  message("\tDimensions of the normalized count table: ", result$nb.genes, " genes x ", result$nb.samples, " samples. ")
  message.with.time("Finished NormalizeCounts() for Recount experiment ID ", parameters$recountID)

  return(result)
}
