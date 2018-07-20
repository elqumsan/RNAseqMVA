#' @title scaling RNA-seq count table
#' @author Jacques van Helden and Mustafa AbuElqumsan
#' @description normalisation of RNA-seq count table.
#' More precisely this function runs a sample-wise scaling so that all the samples
#' have the same value for a user-defined scaling parameter.
#' By default, we use the quantile 0.75 as scaling factor.
#'
#' @param self an object belonging to the class DataTableWithClasses
#' This input object contains a count table, with one column per sample and
#' one row per gene + a phenoTable + sample class specifications.
#'
#' The phenoTable of the input object may be filtered out adapted it in case samples would be suppressed because they have a scaling factor of 0.
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
#' @return the function returns an object of class DataTableWithClasses, where the dataTable
#' coontains the normalised counts. Note that the normalized count table may have a smaller
#' number of samples than the input count table because if some samples have a scaling
#' factor of 0 they will be filtered out. In such case, the phenoTable and classLabels will
#' be filtered out in the same way to ensure consistency  of the attributes of the
#' returned DataTableWithClass object.
#'
#'
#' @examples
#'
#' @export
NormalizeSamples <- function(self,
                             parameters,
                             # method = parameters$standardization$method,
                             # quantile = parameters$standardization$quantile,
                             log2 = TRUE,
                             epsilon = 0.1) {

  if (log2) {
    message.with.time("Starting log2 NormalizeSamples() for Recount experiment ID ", self[["ID"]])
  } else{
    message.with.time("Starting NormalizeSamples() for Recount experiment ID ", self[["ID"]])
  }

  ## Check required parameters
  for (p in c("standardization")) {
    if (is.null(parameters[[p]])) {
      stop("NormalizeSamples()\tMissing required parameter: '", p,
           "'.\n\tPlease check configuration file. ")
    } else {
      assign(p, parameters[[p]])
    }
  }

  if (is.null(parameters$standardization$method)) {
    stop("NormalizeSamples()\tMissing required parameter: standardization method")
  }
  method = parameters$standardization$method

  # method = parameters$standardization$method,
  # quantile = parameters$standardization$quantile,

  # counts <- self$dataTable
  # phenoTable <- self$phenoTable
  # classLabels  <- self$classLabels


  #dim(self$dataTable)
  ## Compute sample-wise statistics
  message("\t", "Computing sample-wise statistics\t", recountID)
  sampleStats <- data.frame(
    sum = apply(self$dataTable, 2, sum, na.rm = TRUE),
    min = apply(self$dataTable, 2, min, na.rm = TRUE),
    mean = apply(self$dataTable, 2, mean, na.rm = TRUE),
    Q1 = apply(self$dataTable, 2, quantile, na.rm = TRUE, probs = 0.25),
    median = apply(self$dataTable, 2, median, na.rm = TRUE),
    Q3 = apply(self$dataTable, 2, quantile, na.rm = TRUE, probs = 0.75),
    max = apply(self$dataTable, 2, max, na.rm = TRUE),
    zero.values = apply(self$dataTable == 0, 2, sum),
    na.values = apply(is.na(self$dataTable), 2, sum),
   # infinite.values = apply(is.infinite(self$dataTable), 2, sum)
   infinite.values = length(is.infinite(unlist(self$dataTable)) == TRUE)
  )

  ## check it the original data contains NA values
  if (sum(sampleStats$na.values) > 0) {
    message("\tRaw counts contain ", sum(sampleStats$na.values), " NA values. ")
  }

  if (method == "quantile") {
    if (is.null(parameters$standardization$quantile)) {
      stop("NormalizeSamples()\tMissing required parameter: standardization quantile")
    }
    quantile = parameters$standardization$quantile

    message("\tNormalizing counts. Method = edgeR RLE.")
    d <- DGEList(counts = self$dataTable, group = self$classLabels)
    # d$samples$group <- relevel(d$samples$group) ## Ensure that condition 2 is considered as the reference
    d <- calcNormFactors(d, method = "upperquartile", p = quantile)                 ## Compute normalizing factors
    sampleStats$scaling.factor <- d$samples$norm.factors
    # message("\tNormalizing counts. Scaling factor: sample quantile ", quantile)
    # sampleStats$norm.quantile <- apply(self$dataTable, 2, quantile, na.rm = TRUE, probs = quantile)
    # sampleStats$scaling.factor <- sampleStats$norm.quantile

  } else if (method == "mean") {
    message("\tNormalizing counts. Scaling factor: sample mean.  ")
    ## Note: mean is very sensitive to outliers, which are very problematic with RNAseq data
  #  sampleStats$scaling.factor <- apply(self$dataTable, 2, quantile, na.rm = TRUE, probs=quantile)
    sampleStats$scaling.factor <- sampleStats$mean

  } else if (method == "median") {
    message("\tNormalizing counts. Scaling factor: sample median.  ")
    sampleStats$scaling.factor <- sampleStats$median

  } else if (method == "sum") {
    message("\tNormalizing counts. Scaling factor: sample sum (= libsum = total counts)  ")
    sampleStats$scaling.factor <- sampleStats$sum

  } else if (method == "TMM") {
    message("\tNormalizing counts. Method = edgeR TMM.")
    d <- DGEList(counts = self$dataTable, group = self$classLabels)
    # d$samples$group <- relevel(d$samples$group) ## Ensure that condition 2 is considered as the reference
    d <- calcNormFactors(d, method = "TMM")                 ## Compute normalizing factors
    sampleStats$scaling.factor <- d$samples$norm.factors

  } else if (method == "RLE") {
    message("\tNormalizing counts. Method = edgeR RLE.")
    d <- DGEList(counts = self$dataTable, group = self$classLabels)
    # d$samples$group <- relevel(d$samples$group) ## Ensure that condition 2 is considered as the reference
    d <- calcNormFactors(d, method = "RLE")                 ## Compute normalizing factors
    sampleStats$scaling.factor <- d$samples$norm.factors

  } else if (method == "DESeq2") {
    message("\tNormalizing counts. Method = DESeq2.")
    stop("NormalizeSamples method DESeq2 NEEDS TO BE IMPLEMENTED")

  } else {
    stop(method, " is not a valid method for NormalizeSamples()")
  }

  ## Detect problems related to null scaling factors, which may happen in some datasets due to a very large number of zeros.
  zeroScaledSamples <- sampleStats$scaling.factor == 0
  discaredSampleNames <- colnames(self$dataTable[, zeroScaledSamples])

  # apply(self$dataTable[, zeroScaledSamples]==0, 2, sum)
  message("\tDiscarding ", sum(zeroScaledSamples), " samples because their scaling factor is null. ")
  message("\tDiscarded samples: ", paste(collapse = "; ", discaredSampleNames))

  ## Compute normalised counts
  normTarget <- median(sampleStats$scaling.factor[!zeroScaledSamples])
  normCounts <- t(t(self$dataTable[, !zeroScaledSamples]) / sampleStats$scaling.factor[!zeroScaledSamples]) * normTarget
  # Check that all normalised coutns have the same scaling factor
  # apply(normCounts, 2, quantile, probs=quantile)
  # dim(self$dataTable)
  # dim(self$phenoTable)
  # length(self$classLabels)
  # dim(normCounts)

  ## Run log2 transformation if required
  if (log2) {
    ## TO BE CHANGED
    ## In principle there is no reason for the log2 transformation to generate NA values.
    ## We should certainly not use na.omit, because this removes any sample that would have any NA value.
    ## Instead, we need to understand why there are NA values in the original count table.
    normCounts <- log2(normCounts + epsilon)
  }

  result <- self
  result$dataTable <- normCounts
  if (sum(zeroScaledSamples) > 0) {
    result$phenoTable <- self$phenoTable[!zeroScaledSamples,]

    ## Update  sample-related parameters
    result$sampleNames <- self$sampleNames[!zeroScaledSamples]
    result$nbSamples <- length(result$sampleNames)

    ## Re-build class-related parameters
    result <- buildAttributes(result)

    ## TO DO !!!!
    ## AFTER HAVING DISCARDD SAMPLES, WE NEED TO CHECK THAT THE NB OF SAMPLES PER CLAS IS STILL ABOVE THE THRESHOLD
  }

  # hist(unlist(normCounts), breaks=100)

  ## Build the result list
 # result <- list()
  result$method <- method
  if (method == "quantile") {
    result$quantile <- quantile
  }
  result$sampleStats <- sampleStats
  result$discaredSampleNames <- discaredSampleNames

  if (log2) {
    result$dataType <- "log2norm_counts"
    message("\tDimensions of the log2 normalized count table: ", result$nbGenes, " genes x ", result$nbSamples, " samples. ")
    message.with.time("Finished log2 NormalizeSamples() for Recount experiment ID ", result[["ID"]])
  } else{
    result$dataType <- "norm_counts"
    message("\tDimensions of the normalized count table: ", result$nbGenes, " genes x ", result$nbSamples, " samples. ")
    message.with.time("Finished NormalizeSamples() for Recount experiment ID ", result[["ID"]])
  }
  return(result)
}
