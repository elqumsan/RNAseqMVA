#' @title Apply various filters on the observations and variables of
#' a count table in order to prepare it for classification.
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @description Apply various filters on the observations (e.g. biological samples)
#' and variables (e.g. genes) of
#' a count table in order to prepare it for classification.
#' @param countTable a data.frame with one row per observation and one column per variable
#' @param phenoTable a data.frame with one row per observation and one column per attribute
#' @param classColumn="tissue" name of a column of the pheno table which contains the class labels (default: "tissue").
#' In some cases classes must be built by concatenating several columns of the pheno table (e.g. "tissue" and "cell.type" for dataset SRP057196),
#' This can be achieved by providing a vector of column names from the pheno table. In this case, class names
## are built by concantenating the values in the specified columns (separated by "_").
#' @param minSamplesPerClass=parameters$minSamplesPerClass minimum nuimber of samples per class to keep
#' @param nearZeroVarFilter=parameters$nearZeroVarFilter if TRUE, applyb caret::nearZeroVariance() to filter out poor predictor genes
#' @param draw.plot=FALSE if TRUE, draw an histogram of variance per gene
#'
#' @examples
#' ## Load a data set
#' countdata <- loadRecountExperiment(recountID = "SRP048759", mergeRuns = TRUE)
#' countTable <- t(countdata$merged$sampleCounts)
#' phenoTable <- countdata$merged$samplePheno
#'
#' # Run the filtering
#' filteredData <- filterCountTable(countTable, phenoTable, classColumn="tissue")
#'
#' # Replace unfiltered data by filtered data
#' countTable <- filteredData$countTable
#' phenoTable <- filteredData$phenoTable
#' classes <- filteredData$classes
#'
#' ## Filter a dataset and build classes based on 2 columns
#' countdata <- loadRecountExperiment(recountID = "SRP057196", mergeRuns = TRUE)
#' countTable <- t(countdata$merged$sampleCounts)
#' phenoTable <- countdata$merged$samplePheno
#'
#' #' # Run the filtering
#' filteredData <- filterCountTable(countTable, phenoTable, classColumn=c("tissue", "cell.type"), minSamplesPerClass=5)
#' table(filteredData$classes)
#'
#' @import recount
#' @import SummarizedExperiment
#' @import caret
#' @import lattice
#' @import ggplot2
#' @export
filterCountTable <- function(countTable,
                             phenoTable,
                             classColumn = parameters$classColumn,
                             minSamplesPerClass = parameters$minSamplesPerClass,
                             nearZeroVarFilter = parameters$nearZeroVarFilter,
                             draw.plot = TRUE) {

  message.with.time("Filtering count table")
  result <- list()

  ## Check if there are NA values, and discard all genes having at least one NA value
  if (sum(is.na(countTable)) > 0) {
    naGenes <- colnames(countTable)[apply(is.na(countTable), 2, sum) > 0]
    result$naGenes <- naGenes
    message("\tFiltering out ", length(naGenes)," genes with NA values")
    if (naGenes > 0) {
      keptGenes <- setdiff(colnames(countTable), naGenes)
      countTable <- countTable[, keptGenes]
      dim(countTable)
    }
    ## Report the number of genes with NA values
    if (exists("naGenes")){
      message("\tThis count Table contains ", length(result$naGenes), " genes with NA values")
      #    result$naGenes
    }

  } else {
    message("\tThis count Table does not contain any NA values. ")
    result$naGenes <- vector()
  }

  ## Detect genes with zero variance
  message("\tDetecting genes with zero variance")
  varPerGene <- apply(countTable, 2 , var, na.rm=TRUE) # Compute variance per gene (row)
  # table(varPerGene == 0)

  # sum(zeroVarGenes == 0) # count genes with zero variance
  zeroVarGenes <- names(varPerGene)[varPerGene == 0] # identify genes having zero variance
  # length(zeroVarGenes) ## Number of genes to discard
  keptGenes <- names(varPerGene)[varPerGene > 0]
  #length(keptGenes) ## Number of genes to keep
  countTable <- countTable[, keptGenes]
  #dim(countTable)
  message("\tFiltering out ", length(zeroVarGenes)," genes with zero variance; keeping ", length(keptGenes), " genes")



  ## Use caret::nearZeroVar() to discard supposedly bad predictor genes.
  ## This includes zero variance genes (already filtered above) but also less trivial cases, see nearZeroVar() doc for details.
  message("\tDetecting genes with near-zero variance  with caret::nearZeroVar()")
  nearZeroVarColumns <- nearZeroVar(countTable, allowParallel = T, saveMetrics = FALSE)
  nearZeroVarGenes <- colnames(countTable)[nearZeroVarColumns]
  #length(nearZeroVarGenes)
  keptGenes <- setdiff(colnames(countTable), nearZeroVarGenes)
  message("\tFiltering out ", length(nearZeroVarGenes)," genes with near-zero variance (poor predictors); kept genes: ", length(keptGenes))
  countTable <- countTable[, keptGenes]
  # dim(countTable)


  ## Count the number of zero values per gene
  # zerosPerGene <- apply(countTable == 0, 2, sum)
  # hist(zerosPerGene)
  # hist(zerosPerGene[nearZeroVarGenes])


  ## Plot an histogram to compare variance distribution  between all genes and those with  near-zero variance
  if (draw.plot) {
    plot.file <- file.path(dir.NormImpact, "var_per_gene_hist.pdf")
    message("\tVariance per gene histograms\t", plot.file)
    pdf(plot.file, width=7, height=8)

    logVarPerGene <- log2(varPerGene)
    noInfVar <- !is.infinite(logVarPerGene)
    xmin <- floor(min(logVarPerGene[noInfVar]))
    xmax <- ceiling(max(logVarPerGene[noInfVar]))
    xlim <- c(xmax, xmin)
    breaks <- seq(from=xmin,  to=xmax, by=0.1)
    par(mfrow=c(3,1))
    hist(log2(varPerGene[noInfVar]),
         breaks=breaks,
         col="gray", border = "gray",
         main = paste("All non-zero var genes; ", parameters$recountID),
         xlab="log2(varPerGene)",
         ylab="Number of genes")
#    legend("topright", parameters$recountID)
    hist(log2(varPerGene[nearZeroVarGenes]),
         breaks=breaks,
         col="red", border = "red",
         main = "Near zero variance",
         xlab="log2(varPerGene)",
         ylab="Number of genes")
    hist(log2(varPerGene[keptGenes]),
         breaks=breaks,
         col="darkgreen", border = "darkgreen",
         main = "Kept genes",
         xlab="log2(varPerGene)",
         ylab="Number of genes")
    par(mfrow=c(1,1))

    silence <- dev.off()
  }

  # sample.nb <- nrow(countTable)
  # gene.nb <- ncol(countTable)
  # dim(countTable)
  # dim(phenoTable)
  # View( phenoTable)

  ################################################################
  ## Specify sample classes by extracting information about specified class columns
  if (is.null(classColumn) || (length(classColumn) < 1)) {
    stop("classColumn must be defined. ")
  } else if (length(classColumn) == 1) {
    classes <-  as.vector(phenoTable[, classColumn])
  } else {
    ## Combine several columns to establish the classes
    classes <- apply(phenoTable[, classColumn], 1, paste, collapse="_")
  }
  # table(classes)


  ################################################################
  ## Check if there are NA values in the sample classes
  discardedSamples <- is.na(classes)
  if (sum(discardedSamples) > 0) {
    message("\tDiscarding ", sum(discardedSamples), " samples with undefined class in ", classColumn, " column of pheno table")
    nonaSamples <- !is.na(tissue)
    countTable <- countTable[nonaSamples, ]
    phenoTable<- phenoTable[nonaSamples, ]
    classes <- classes[nonaSamples]
    sample.nb <- nrow(countTable)
  }

  ################################################################
  ## Select classess for which we dispose of at least 10 samples
  # length(classes)
  message("\tSelecting classes with at least ", minSamplesPerClass, " samples")
  samplesPerClass <- table(classes)
  discardedClasses <- names(samplesPerClass)[samplesPerClass < minSamplesPerClass]
  selectedClasses <- names(samplesPerClass)[samplesPerClass >= minSamplesPerClass]
  selectedSamples <- classes %in% selectedClasses
  if (length (discardedClasses) > 0) {
    message("\tDiscarding ", length (discardedClasses), " classes containng less than ", minSamplesPerClass, "samples")
    message("\tDiscarded classes\t", paste(collapse=",", discardedClasses))
  }
  message("\tKeeping ", sum(selectedSamples), " samples from ",
          length(selectedClasses), " classes")
  message("\tKept classes\t", paste(selectedClasses, collapse=", "))


  ## Update count table, pheno table and classes vector to keep only the samples belonging to selected classess
  countTable <- countTable[selectedSamples, ]
  phenoTable<- phenoTable[selectedSamples, ]
  classes <- classes[selectedSamples]
  # sample.nb <- nrow(countTable)

  ################################################################################
  ## Return an object with the filtered counts + updated pheno table selected classes and samples
  ## Transpose the table in order to get it in the suitable format for classifiers:
  ## one row per individual, one column per variable.
  message("\tAfter filtering, count table contains ",
          nrow(countTable), " samples (rows) and ",
          ncol(countTable), " genes (columns) ",
          "belonging to ", length(selectedClasses), " classes")
  # View(countTable)
  # dim(countTable)

  result$countTable <- countTable
  result$phenoTable <- phenoTable
  result$classes <- classes
  result$selectedClasses <- selectedClasses

  result$zeroVarGenes <- zeroVarGenes
  result$nearZeroVarGenes <- nearZeroVarGenes
  result$keptGenes <- keptGenes
  message.with.time("Finished Filter Count Table process for Recount experiment ID ", parameters$recountID)

  return(result)

}
