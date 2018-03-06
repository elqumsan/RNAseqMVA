#' @title Apply various filters on the observations and variables of
#' a count table in order to prepare it for classification.
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @description Apply various filters on the observations (e.g. biological samples)
#' and variables (e.g. genes) of
#' a count table in order to prepare it for classification.
#' @param countTable a data.frame with one row per observation and one column per variable
#' @param phenoTable a data.frame with one row per observation and one column per attribute
#' @param classColumn=parameters$classColumn name of a column of the pheno table which contains the class labels.
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
  filteredCountTable <- countTable[, keptGenes]
  # dim(countTable)
  # dim(filteredCountTable)
  message("\tFiltering out ", length(zeroVarGenes)," genes with zero variance; keeping ", length(keptGenes), " genes")



  ## Use caret::nearZeroVar() to discard supposedly bad predictor genes.
  ## This includes zero variance genes (already filtered above) but also less trivial cases, see nearZeroVar() doc for details.
  if (nearZeroVarFilter) {

    message("\tDetecting genes with near-zero variance using caret::nearZeroVar()")
    nearZeroVarColumns <- nearZeroVar(countTable, allowParallel = T, saveMetrics = FALSE)
    nearZeroVarGenes <- colnames(countTable)[nearZeroVarColumns]
    #length(nearZeroVarGenes)
    keptGenes <- setdiff(colnames(filteredCountTable), nearZeroVarGenes)
    message("\tFiltering out ", length(nearZeroVarGenes)," genes with near-zero variance (poor predictors); kept genes: ", length(keptGenes))
    filteredCountTable <- filteredCountTable[, keptGenes]
  # dim(filteredCountTable)
  }


  ## Plot an histogram to compare variance distribution  between all genes and those with  near-zero variance
  if (draw.plot) {
    plot.file <- file.path(dir.NormImpact, "var_per_gene_hist.pdf")
    message("\tVariance per gene histograms\t", plot.file)
    pdf(plot.file, width=7, height=12)

    logVarPerGene <- log2(varPerGene)
    noInfVar <- !is.infinite(logVarPerGene)
    xmin <- floor(min(logVarPerGene[noInfVar]))
    xmax <- ceiling(max(logVarPerGene[noInfVar]))
    xlim <- c(xmax, xmin)
    varbreaks <- seq(from=xmin,  to=xmax, by=0.1)
    if (nearZeroVarFilter) {
      par(mfrow=c(4,1))
    } else {
      par(mfrow=c(3,1))
    }
    hist(log2(varPerGene[noInfVar]),
         breaks=varbreaks,
         col="gray", border = "gray",
         main = paste("All non-zero var genes; ", parameters$recountID),
         xlab="log2(varPerGene)",
         ylab="Number of genes")
#    legend("topright", parameters$recountID)
    if (nearZeroVarFilter) {
      hist(log2(varPerGene[nearZeroVarGenes]),
         breaks=varbreaks,
         col="red", border = "orange",
         main = "Near zero variance",
         xlab="log2(varPerGene)",
         ylab="Number of genes")
    }
    hist(log2(varPerGene[keptGenes]),
         breaks=varbreaks,
         col="darkgreen", border = "#00BB00",
         main = "Kept genes",
         xlab="log2(varPerGene)",
         ylab="Number of genes")

    ## Count the number of zero values per gene
    zerosPerGene <- apply(countTable == 0, 2, sum)
    zerobreaks <- seq(from=0, to=max(zerosPerGene+1), length.out =50)
    # zerobreaks <- seq(from=0, to=max(zerosPerGene+1), by=1)

    #### Histogram of zero values per gene. ####
    ##
    ## Displayed in green (color for kept genes) because we will then
    ## overlay the histograms of near zero var (orange) and zero var (red),
    ## so that the remaining part of the histogram will correspond to genes
    ## kept.
    hist(zerosPerGene,
         breaks = zerobreaks,
         main = "Zeros per gene",
         xlab = "Number of zero values",
         ylab = "Number of genes",
         col = "#00BB00", border = "#00BB00")
    if (nearZeroVarFilter) {
      hist(zerosPerGene[union(nearZeroVarGenes, zeroVarGenes)],
         breaks = zerobreaks,
         add=TRUE, col="orange", border="orange")
    }
    hist(zerosPerGene[zeroVarGenes],
         breaks = zerobreaks,
         add=TRUE, col="red", border="red")
    if (nearZeroVarFilter) {
      legend("topleft",
             legend=paste(
               sep= ": ",
               c("Kept genes", "Near-zero variance", "Zero variance"),
               c(length(keptGenes), length(nearZeroVarGenes), length(zeroVarGenes))),
             lwd=5,
             col=c("#00BB00", "orange", "red")
      )
    } else {
      legend("topleft",
             legend=paste(
               sep= ": ",
               c("Kept genes",  "Zero variance"),
               c(length(keptGenes), length(zeroVarGenes))),
             lwd=5,
             col=c("#00BB00", "red"))

    }
    #    hist(zerosPerGene[nearZeroVarGenes])

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
    nonaSamples <- !is.na(classes)
    filteredCountTable <- filteredCountTable[nonaSamples, ]
    # dim(filteredCountTable)
    # dim(filteredPhenoTable)
    filteredPhenoTable<- phenoTable[nonaSamples, ]
    filteredClasses <- classes[nonaSamples]
  } else {
    filteredPhenoTable<- phenoTable
    filteredClasses <- classes
  }
  filterdSampleNb <- nrow(filteredCountTable)

  ################################################################
  ## Select classess for which we dispose of at least 10 samples
  # length(classes)
  message("\tSelecting classes with at least ", minSamplesPerClass, " samples")
  samplesPerClass <- table(filteredClasses)
  discardedClasses <- names(samplesPerClass)[samplesPerClass < minSamplesPerClass]
  selectedClasses <- names(samplesPerClass)[samplesPerClass >= minSamplesPerClass]
  selectedSamples <- filteredClasses %in% selectedClasses
  if (length (discardedClasses) > 0) {
    message("\tDiscarding ", length (discardedClasses), " classes containng less than ", minSamplesPerClass, "samples")
    message("\tDiscarded classes\t", paste(collapse=",", discardedClasses))
  }
  message("\tKeeping ", sum(selectedSamples), " samples from ",
          length(selectedClasses), " classes")
  message("\tKept classes\t", paste(selectedClasses, collapse=", "))


  ## Update count table, pheno table and classes vector to keep only the samples belonging to selected classess
  filteredCountTable <- filteredCountTable[selectedSamples, ]
  # dim(filteredCountTable)
  filteredPhenoTable <- filteredPhenoTable[selectedSamples, ]
  # dim(filteredPhenoTable)
  filteredClasses <- filteredClasses[selectedSamples]
  # length(filteredClasses)
  filteredSampleNb <- nrow(filteredCountTable)

  ################################################################################
  ## Return an object with the filtered counts + updated pheno table selected classes and samples
  ## Transpose the table in order to get it in the suitable format for classifiers:
  ## one row per individual, one column per variable.
  message("\tBefore filtering, count table contains ",
          nrow(countTable), " samples (rows) and ",
          ncol(countTable), " genes (columns) ",
          "belonging to ", length(classes), " classes")
  message("\tAfter filtering, count table contains ",
          nrow(filteredCountTable), " samples (rows) and ",
          ncol(filteredCountTable), " genes (columns) ",
          "belonging to ", length(selectedClasses), " classes")
  # View(filteredCountTable)
  # dim(filteredCountTable)


  ## Return unfiltered count table + phenotable + classes
  result$countTable <- countTable
  result$phenoTable <- phenoTable
  result$classes <- classes

  ## Include the filtering criteria in the returned list
  result$zeroVarGenes <- zeroVarGenes
  if (nearZeroVarFilter) {
    result$nearZeroVarGenes <- nearZeroVarGenes
  }
  result$keptGenes <- keptGenes
  result$selectedClasses <- selectedClasses

  ## Return gene-wise and sample-wise filtered count table
  result$filteredCountTable <- filteredCountTable
  result$filteredPhenoTable <- filteredPhenoTable
  result$filteredClasses <- filteredClasses

  message.with.time("Finished filterCountTable() for Recount experiment ID ", parameters$recountID)

  return(result)

}
