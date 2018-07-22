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
#' @param method="quantile" normalization method.
#' Supported normalization methods: sum, mean, median, quantile, TMM, RLE, DESeq2.
#'
#' Sum and mean are really not recommended because very sensitive to outliers.
#' They are implemented only for the sake of comparison.
#'
#' Quantile: Preferred quantile = 0.75. Quantiles are more robust to outliers than the mean
#' or sum. The median is sometimes weak in RNA-seq data, the 75th percentile (quantile 0.75)
#' seems a good tradeoff between robustness and representativity (taking into account a
#' representative proportion of the total counts per samples).
#'
#' Median: actually runs quantile-based scaling factor with quantile=0.5.
#'
#' TMM: trimmed mean of M-values proposed by Robinson and Oshlack (2010), computed via edgeR::calcNormFactors().
#'
#' RLE: relative log expression proposed by Anders and Huber (2010), computed via edgeR::calcNormFactors().
#'
#' DESeq2: compute size factors via DESeq2::estimateSizeFactors()
#'
#' @param quantile=0.75 quantile used as scaling factor when quantile method is selected.
#'
#' @param log2=FALSE apply a log2 transformation
#'
#' @param epsilon=0.1 value added to all counts before applying the log2 transformation
#' in order to avoid zero counts.
#'
#' @param detailed.sample.stats=FALSE compute detailed sample stats (takes some seconds)
#'
#' @return the function returns an object of class DataTableWithClasses, where the dataTable
#' contains the normalised counts. Note that the normalized count table may have a smaller
#' number of samples than the input count table because if some samples have a scaling
#' factor of 0 they will be filtered out. In such case, the phenoTable and classLabels will
#' be filtered out in the same way to ensure consistency  of the attributes of the
#' returned DataTableWithClass object.
#'
#'
#' @examples
#' @export
NormalizeSamples <- function(self,
                             method = "quantile",
                             quantile = 0.75,
                             log2 = FALSE,
                             epsilon = 0.1,
                             detailed.sample.stats = FALSE) {

  ## Start message
  # if (log2) {
  #   message.with.time("Starting log2 NormalizeSamples() for Recount experiment ID ", self[["ID"]])
  # } else{
  #   message.with.time("Starting NormalizeSamples() for Recount experiment ID ", self[["ID"]])
  # }

  ## Define method name
  if (method == "quantile") {
    if (!exists("quantile")) {
      stop("NormalizeSamples()\tMissing required parameter: standardization quantile")
    }
    if (is.null(quantile)) {
      stop("NormalizeSamples()\tquantile-based scaling requires a non-null quantile parameter.")
    }
    method.name <- paste(sep = "", "q", quantile)
  } else {
    method.name <- method
  }

  ## Median will be treated as quantile 0.5
  if (method == "median") {
    quantile <- 0.5
  }

  ## Append log2 suffix if log2-transformation is activated
  if (log2) {
    method.name <- paste(sep = "", method.name, "_log2")
  }
  message("\t", "Normalizing\t", recountID, "\tmethod: ", method.name)

  ## Extract the count matrix
  counts <- self$dataTable
  parameters <- self$parameters

  ## if required, discard the zero counts before computing size factors
  if (is.null(parameters$standardization$nozero)) {
    ## By default, nozero parameter is activated
    nozero <- TRUE
    # View(head(non.null.counts))
  } else {
    nozero <- parameters$standardization$nozero
  }
  counts.nozero <- counts
  counts.nozero[counts  == 0] <- NA
  if (nozero) {
    message("\t\tIgnoring zero values for normalization")
    # View(head(non.null.counts))
    counts <- counts.nozero
  } else {
    counts <- self$dataTable
  }


  ## Compute sample-wise statistics
  message("\t\t", "Computing sample-wise statistics\t", recountID)
  sampleStats <- data.frame(
    zero.values = apply(self$dataTable == 0, 2, sum),
    na.values = apply(is.na(self$dataTable), 2, sum),
    sum = apply(self$dataTable, 2, sum, na.rm = TRUE),
    mean = apply(self$dataTable, 2, mean, na.rm = TRUE),
    mean.nozero = apply(counts.nozero, 2, mean, na.rm = TRUE))
  # head(sampleStats)
  ## Report number of NA and zero values
  message("\t\tOrginal data table contains ",
          format(x = sum(sampleStats$zero.values), big.mark = ","),
          " zeros ")
  if (sum(sampleStats$na.values) > 0) {
    message("\t\tRaw counts contain ",
            format(x = sum(sampleStats$na.values), big.mark = ","),
            " NA values. ")
  }

  if (detailed.sample.stats) {
    sampleStats$min <- apply(self$dataTable, 2, min, na.rm = TRUE)
    sampleStats$Q1 <- apply(self$dataTable, 2, quantile, na.rm = TRUE, probs = 0.25)
    sampleStats$median <- apply(self$dataTable, 2, median, na.rm = TRUE)
    sampleStats$Q3 <- apply(self$dataTable, 2, quantile, na.rm = TRUE, probs = 0.75)
    sampleStats$max <- apply(self$dataTable, 2, max, na.rm = TRUE)
    sampleStats$infinite.values <- apply(is.infinite(as.matrix(self$dataTable)), 2, sum)

    ## Compute statistics of zero-ommited data table
    if (nozero) {
      sampleStats.nozero <- data.frame(
        sum = apply(counts, 2, sum, na.rm = TRUE),
        mean = apply(counts, 2, mean, na.rm = TRUE),
        min = apply(counts, 2, min, na.rm = TRUE),
        Q1 = apply(counts, 2, quantile, na.rm = TRUE, probs = 0.25),
        median = apply(counts, 2, median, na.rm = TRUE),
        Q3 = apply(counts, 2, quantile, na.rm = TRUE, probs = 0.75),
        max = apply(counts, 2, max, na.rm = TRUE),
        zero.values = apply(counts == 0, 2, sum),
        na.values = apply(is.na(counts), 2, sum),
        infinite.values = apply(is.infinite(as.matrix(counts)), 2, sum)
      )
      self$sampleStats.nozero <- sampleStats.nozero

      ## Report number of NA and zero values
      message("\t\tZero-ommited data-table contains ",
              format(x = sum(sampleStats.nozero$na.values), big.mark = ","),
              " NA values. ")
    }
  }
  self$sampleStats <- sampleStats
  # head(self$sampleStats)


  if (method %in% c("quantile", "median")) {

    quantile.method <- "RNAseqMVA"
    if (quantile.method == "edgeR") {

      ## Compute quantile-based scaling factor via edgeR
      ## NOTE (2018-07-21) : with single-cell data containing MANY zeros, this returns Inf scaling factors for almost all the samples
      message("\t\tNormalizing counts with edgeR::calcNormFactors(method=upperquartile, p=", quantile,")")
      d <- DGEList(counts = self$dataTable, group = self$classLabels)
      # d$samples$group <- relevel(d$samples$group)
      d <- calcNormFactors(d, method = "upperquartile", p = quantile)                 ## Compute normalizing factors
      sampleStats$scaling.factor <- d$samples$norm.factors
      sampleStats$size.factor <- 1/sampleStats$scaling.factor

    } else {
      message("\t\tScaling factor: sample quantile ", quantile)
      sampleStats$norm.quantile <- apply(counts, 2, quantile, na.rm = TRUE, probs = quantile)
      sampleStats$size.factor <-  sampleStats$norm.quantile
      sampleStats$scaling.factor <- 1 / sampleStats$norm.quantile
      sampleStats$scaling.factor <- sampleStats$scaling.factor / mean(sampleStats$scaling.factor[!is.infinite(sampleStats$scaling.factor)])
      # mean(sampleStats$scaling.factor[!is.infinite(sampleStats$scaling.factor)])
      # hist(sampleStats$scaling.factor, breaks = 1000)
    }
    null.scaling <- sum(sampleStats$size.factor == 0)
    # inf.scaling <- sum(is.infinite(sampleStats$scaling.factor))
    if (null.scaling > 1) {
      message("\t\tdiscarding ", null.scaling, " samples with null value for quantile ", quantile)
    }

  } else if (method %in% c("mean", "libsum", "TC", "sum")) {
    message("\t\tScaling factor: library size (equivalent for sum, mean, total counts).  ")
    ## Note: mean is very sensitive to outliers, which are very problematic with RNAseq data
  #  sampleStats$scaling.factor <- apply(counts, 2, quantile, na.rm = TRUE, probs=quantile)
    sampleStats$size.factor <- sampleStats$mean
    sampleStats$scaling.factor <- 1 / sampleStats$size.factor
    sampleStats$scaling.factor <- sampleStats$scaling.factor / mean(sampleStats$scaling.factor)
    # hist(sampleStats$scaling.factor, breaks = 1000)
    # mean(sampleStats$scaling.factor)


  } else if (method == "TMM") {
    message("\t\tRunning edgeR::calcNormFactors(method='TMM').")
    d <- DGEList(counts = self$dataTable, group = self$classLabels)
    # d$samples$group <- relevel(d$samples$group) ## Ensure that condition 2 is considered as the reference
    d <- calcNormFactors(d, method = "TMM")                 ## Compute normalizing factors
    sampleStats$scaling.factor <- d$samples$norm.factors
    sampleStats$size.factor <- 1/sampleStats$scaling.factor

  } else if (method == "RLE") {
    ## Run edgeR to compute the relative log expression defined by Anders and Huber (2010).
    message("\t\tRunning edgeR::calcNormFactors(method='RLE').")
    d <- DGEList(counts = self$dataTable, group = self$classLabels)
    # d$samples$group <- relevel(d$samples$group) ## Ensure that condition 2 is considered as the reference
    d <- calcNormFactors(d, method = "RLE")                 ## Compute normalizing factors
    sampleStats$scaling.factor <- d$samples$norm.factors
    sampleStats$size.factor <- 1/sampleStats$scaling.factor

  } else if (method == "DESeq2") {
    message("\t\tRunning DESeq2::estimateSizeFactors()")
    ## Replace non-alphanumeric charactersi by "_" in class labels, to avoid message from DESeq2
    classLabels <- gsub(pattern = "[^[:alnum:] ]", replacement = "_", self$classLabels)
    classLabels <- gsub(pattern = " ", replacement = "_", classLabels)
    unique(classLabels)
    dds <- DESeqDataSetFromMatrix(
      countData = self$dataTable,
      colData = data.frame(classes = classLabels), ~ classes)
    dds <- estimateSizeFactors(dds)
    sampleStats$size.factor <- sizeFactors(dds)
    sampleStats$scaling.factor <- 1/sampleStats$size.factor
    # plot(sampleStats$size.factor, sampleStats$sum) ## THE DIFFERENCE IS QUITE IMPRESSIVE

  } else if (method == "VSD") {
    # Compute variance stabilizing transformations (VST) via DESeq2 (Tibshirani 1988; Huber et al. 2003; Anders and Huber 2010)
    message("\t\tRunning DESeq2::vsd()")
    stop("NOT FINISHED THE IMPLEMENTATION YET")
    ## Create a DESeqDataset object from the count table
    dds <- DESeqDataSetFromMatrix(counts = self$dataTable, colData = data.frame(classes = self$classLabels), ~ classes  )
    vsd <- vst(dds, blind = FALSE)
    names(vsd)
    ## QUESTION: HOW DO I GET THE SIZE FACTORS FROM THE RESULTING OBJECT ?

#    rld <- rlog(dds, blind = FALSE)

    # View(vsd)
    # head(assay(vsd), 3)
    # library("pheatmap")
    # library("vsn")
    # meanSdPlot(assay(vsd))

    ## Run  differential expression analysis with DESeq2
#    dds <- DESeq(dds)


#    rld <- rlog(dds, blind=FALSE)

  } else {
    stop(method, " is not a valid method for NormalizeSamples()")
  }

  ## Detect problems related to null scaling factors, which may happen in some datasets due to a very large number of zeros.
  discardedSamples <-
    (sampleStats$size.factor == 0) |
    is.infinite(sampleStats$scaling.factor) |
    is.na(sampleStats$size.factor)
  discaredSampleNames <- vector()
  if (sum(discardedSamples) > 0) {
    discaredSampleNames <- colnames(counts[, discardedSamples])
    message("\t\tDiscarding ", sum(discardedSamples), " samples because their size factor is null or NA. ")
    message("\t\tDiscarded samples: ", paste(collapse = "; ", discaredSampleNames))
  }

  ## Compute normalised counts
  normTarget <- mean(sampleStats$scaling.factor[!discardedSamples]) ## Ensure library eize equality before and after standardization
  sampleStats$scaling.factor <- sampleStats$scaling.factor * normTarget
  normCounts <- t(t(self$dataTable[, !discardedSamples]) * sampleStats$scaling.factor)

  ## Run log2 transformation if required
  if (log2) {
    normCounts <- log2(normCounts + epsilon)
  }

  ## Build the result object
  result <- self ## Clone the input object
  result$dataType <- method.name
  result$method <- method
  if (method == "quantile") {
    result$quantile <- quantile
  }
  result$nozero <- nozero
  result$sampleStats <- sampleStats
  result$discaredSampleNames <- discaredSampleNames
  result$dataTable <- normCounts

  ## In case samples have been discarded, re-build the object to
  ## ensure consistency between sample- and class-dependent attributes
  if (sum(discardedSamples) > 0) {
    result$phenoTable <- self$phenoTable[!discardedSamples,]

    ## Update  sample-related parameters
    result$sampleNames <- self$sampleNames[!discardedSamples]
    result$nbSamples <- length(result$sampleNames)

    ## Re-build class-related parameters
    result <- buildAttributes(result)
  }

  message("\t\tDimensions before normalization: ", format(big.mark = ",", nrow(self$dataTable)), " features x ", ncol(self$dataTable), " samples. ")
  message("\t\tDimensions after normalization: ", format(big.mark = ",", result$nbGenes), " features x ", result$nbSamples, " samples. ")
  # message.with.time("\tFinished NormalizeSamples() for Recount experiment ID ", result[["ID"]])
  return(result)
}
