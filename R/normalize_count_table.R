#' @title Normalize RNA-seq count table
#' @author Jacques van Helden
#' @description Normalize a table of raw counts from RNA-seq
#' @param rawCounts a raw count table, with one column per sample and one row per gene. 
#' Normalization will be done by samples (columns).
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
NormalizeCounts <- function(rawCounts,
                            method="quantile",
                            quantile=0.75,
                            log2=TRUE,
                            epsilon=0.1) {
  
  ## Compute sample-wise statistics
  sampleStats <- data.frame(
    sum=apply(rawCounts, 2, sum, na.rm=TRUE),
    min=apply(rawCounts, 2, min, na.rm=TRUE),
    mean=apply(rawCounts, 2, mean, na.rm=TRUE),
    Q1=apply(rawCounts, 2, quantile, na.rm=TRUE, probs=0.25),
    median=apply(rawCounts, 2, median, na.rm=TRUE),
    Q3=apply(rawCounts, 2, quantile, na.rm=TRUE, probs=0.75),
    max=apply(rawCounts, 2, max, na.rm=TRUE)
  )
  
  if (method == "quantile") {
    sampleStats$norm.quantile <- apply(rawCounts, 2, quantile, na.rm=TRUE, probs=quantile)
    sampleStats$scaling.factor <- sampleStats$norm.quantile
  } else if (method == "mean") {
  #  sampleStats$scaling.factor <- apply(rawCounts, 2, quantile, na.rm=TRUE, probs=quantile)
    sampleStats$scaling.factor <- sampleStats$mean
  } else if (method == "median") {
    sampleStats$scaling.factor <- sampleStats$median
  } else if (method == "sum") {
    sampleStats$scaling.factor <- sampleStats$sum
  } else {
    stop(method, " is not a valid method for NormalizeCounts()")
  }
  
  ## Compute normalised counts
  normTarget <- median(sampleStats$scaling.factor)
  normCounts <- t(t(rawCounts) / sampleStats$scaling.factor) * normTarget
  
  
  ## Run log2 transformation if required
  if (log2) {
    normCounts <- log2(normCounts + epsilon)
  }
  
  # hist(unlist(normCounts), breaks=100)
  
  result <- list()
  result$normCounts <- normCounts
  result$method <- method
  if (method == "quantile") {
    result$quantile <- quantile
  }
  result$sampleStats <- sampleStats
  return(result)
}