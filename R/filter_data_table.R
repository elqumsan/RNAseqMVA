#' @title Apply various filters on the observations and variables of
#' a count table in order to prepare it for classification.
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @description Apply various filters on the individuals (e.g. biological samples)
#' and features (e.g. genes, transcripts, ...) of a count table in order to prepare
#' it for classification.
#' @param countsWithClasses an object of the class DataTableWithClasses.
#' This object contains the data table + the parameters (including filtering parameters).
#' @param draw.plot=TRUE if TRUE, draw an histogram of variance per feature.
#' @param plot.heigh=NULL htigh of the pdf plot. If NULL, computed autopatically depending on the number of panels.
#'
#' @examples
#' ## Load a data set
#' countdata <- loadRecountExperiment(recountID = "SRP048759", mergeRuns = TRUE)
#' dataTable <- t(countdata$merged$sampleCounts)
#' phenoTable <- countdata$merged$samplePheno
#'
#' # Run the filtering
#' filteredData <- filterDataTable(dataTable, phenoTable, classColumn="tissue")
#'
#' # Replace unfiltered data by filtered data
#' dataTable <- filteredData$dataTable
#' phenoTable <- filteredData$phenoTable
#' classLabels <- filteredData$classLabels
#'
#' ## Filter a dataset and build classLabels based on 2 columns
#' countdata <- loadRecountExperiment(recountID = "SRP057196", mergeRuns = TRUE)
#' dataTable <- t(countdata$merged$sampleCounts)
#' phenoTable <- countdata$merged$samplePheno
#'
#' #' # Run the filtering
#' filteredData <- filterDataTable(dataTable, phenoTable, classColumn=c("tissue", "cell.type"), minSamplesPerClass=5)
#' table(filteredData$classLabels)
#'
#' @import lattice
#' @import ggplot2
#' @import caret
#' @export
filterDataTable <- function(rawCounts,
                            draw.plot = TRUE) {

  LoadRequiredBioconductorPackages(c("recount", "SummarizedExperiment"))


  message.with.time("Filtering count table")

  #### Check validity of the input class ####
  if (!is(rawCounts, "DataTableWithClasses")) {
    stop("filterDataTable() requires an object of class DataTableWithClasses")
  }

  #### Check parameters ####


  ## Take parameters from the DataTableWithClasses object
  parameters <- rawCounts$parameters

  ## Check required parameters
  if (is.null(parameters$filtering$na.rm)) {
    message("\tfilterDataTable()\tundefined filtering parameter: na.rm. Setting to default (TRUE).")
    parameters$filtering$na.rm <- TRUE
  }
  na.rm <- parameters$filtering$na.rm

  ## Check min samples per class
  if (is.null(parameters$filtering$minSamplesPerClass)) {
    message("\tfilterDataTable()\tundefined filtering parameter: minSamplesPerClass. Setting to default (10).")
    parameters$filtering$minSamplesPerClass <- 10
  }
  minSamplesPerClass <- parameters$filtering$minSamplesPerClass

  if (is.null(parameters$filtering$nearZeroVarFilter)) {
    message("\tfilterDataTable()\tundefined filtering parameter: nearZeroVarFilter. Setting to default (FALSE).")
    parameters$filtering$nearZeroVarFilter <- FALSE
  }
  nearZeroVarFilter <- parameters$filtering$nearZeroVarFilter

  ## Set the updated parameters to the object
  rawCounts$parameters <- parameters

  ## Initialise a vector with the features that will be kept through the different filtering stpes.
  keptGenes <- rawCounts$featureNames


  #### Treat NA values ####
  if (sum(is.na(rawCounts$dataTable)) > 0) {
    naGenes <- rownames(rawCounts$dataTable)[apply(is.na(rawCounts$dataTable), 1, sum) > 0]
    ## Report the number of genes with NA values
    message("\tThis count Table contains ", length(naGenes), " genes with NA values")
    message("\tFiltering out ", length(naGenes)," genes with NA values")
  } else {
    message("\tThis count Table does not contain any NA values. ")
    naGenes <- vector()
  }

  if (na.rm) {
    keptGenes <- setdiff(rawCounts$featureNames, naGenes)
    message("\tNA filter: discarding ", length(naGenes), " genes out of ", rawCounts$nbGenes, "; keeping ", length(keptGenes))
  }
  #filteredDataTable <- rawCounts$dataTable[, keptGenes]


  #### Detect genes with zero variance ####
  message("\tDetecting genes with zero variance")
  varPerGene <- apply(rawCounts$dataTable, 1 , var, na.rm = TRUE) # Compute variance per gene (row)
  # table(varPerGene == 0)
  zeroVarGenes <- names(varPerGene)[varPerGene == 0] # identify genes having zero variance
  # length(zeroVarGenes) ## Number of genes to discard
  keptGenes <- setdiff(keptGenes, zeroVarGenes)
  # length(keptGenes) ## Number of genes to keep
  message("\tZero var filter: discarding ", length(zeroVarGenes)," genes with zero variance; keeping ", length(keptGenes), " genes")

  #### Use caret::nearZeroVar() to discard supposedly bad predictor genes. #####
  ## This includes zero variance genes (already filtered above) but also less trivial cases, see nearZeroVar() doc for details.
  nearZeroVarGenes <- NULL
  if (nearZeroVarFilter) {
    message("\tDetecting genes with near-zero variance using caret::nearZeroVar()")
    nearZeroVarIndices <- nearZeroVar(t(rawCounts$dataTable), allowParallel = TRUE, saveMetrics = FALSE)
    nearZeroVarGenes <- rownames(rawCounts$dataTable)[nearZeroVarIndices]
    nearZeroVarGenes <- setdiff(nearZeroVarGenes, zeroVarGenes)
    keptGenes <- setdiff(keptGenes, nearZeroVarGenes)
    message("\tNear zero var filter: discarding ", length(nearZeroVarGenes)," genes with near-zero variance (poor predictors); kept genes: ", length(keptGenes))
  }

  #### Extract the table of counts with the subset of kept genes at the end of the different gene filters. ####
  filteredDataTable <- rawCounts$dataTable[keptGenes,]
  # dim(filteredDataTable)

  #### Check if there are NA values in the sample classes ####
  discardedSamples <- is.na(rawCounts$classLabels)
  if (sum(discardedSamples) > 0) {
    message("\tDiscarding ", sum(discardedSamples), " samples with undefined class in ", classColumn, " column of pheno table")
    nonaSamples <- !is.na(rawCounts$classLabels)
    filteredDataTable <- filteredDataTable[, nonaSamples]
    # dim(filteredDataTable)
    # dim(filteredPhenoTable)
    filteredPhenoTable <- rawCounts$phenoTable[nonaSamples, ]
    filteredClassLabels <- rawCounts$classLabels[nonaSamples]
    # table(countsWithClasses$classLabels, filteredClassLabels)
  } else {
    filteredPhenoTable <- rawCounts$phenoTable
    filteredClassLabels <- rawCounts$classLabels
  }

  #### Discard classes having less samples than the user-specified threshold ####
  # length(countsWithClasses$classLabels)
  message("\tSelecting classes with at least ", minSamplesPerClass, " samples")
  samplesPerClass <- table(filteredClassLabels)
  discardedClasses <- names(samplesPerClass)[samplesPerClass < minSamplesPerClass]
  keptClasses <- names(samplesPerClass)[samplesPerClass >= minSamplesPerClass]
  keptSamples <- filteredClassLabels %in% keptClasses
  if (length(discardedClasses) > 0) {
    message("\tDiscarding ", length(discardedClasses), " classes containing less than ", minSamplesPerClass, " samples")
    message("\t\tDiscarded classes\t", paste(collapse = ",", discardedClasses))
  }
  message("\tKeeping ", sum(keptSamples), " samples from ",
          length(keptClasses), " classes")
  message("\tKept classes\t", paste(keptClasses, collapse = ", "))


  ## Update count table, pheno table and countsWithClasses$classLabels vector to keep only the samples belonging to selected classess
  filteredDataTable <- filteredDataTable[, keptSamples]
  # dim(filteredDataTable)
  filteredPhenoTable <- filteredPhenoTable[keptSamples, ]
  # dim(filteredPhenoTable)
  filteredClassLabels <- filteredClassLabels[keptSamples]
  # length(filteredClassLabels)


  ## Return unfiltered count table + phenotable + countsWithClasses$classLabels
  result <- DataTableWithClasses(dataTable = filteredDataTable,
                                 phenoTable =  filteredPhenoTable,
                                 # classColumn = rawCounts$classColumn,
                                 # classColors = rawCounts$classColors,
                                 # variablesType = "all",
                                 dataType = "filtered_counts",
                                 parameters = rawCounts$parameters)

  # class(result)
  # summary(result)


  ## Include the filtering criteria in the returned list
  result$naGenes <- naGenes
  result$zeroVarGenes <- zeroVarGenes
  result$nearZeroVarGenes <- nearZeroVarGenes
  result$keptGenes <- keptGenes
  result$keptClasses <- keptClasses
  result$varPerGeneRaw <- varPerGene
  result$varPerGeneFiltered <- apply(result$dataTable, 1 , var, na.rm = TRUE) # Compute variance per gene (row)

  ## Plot an histogram to compare variance distribution  between all genes and those with  near-zero variance
  if (draw.plot) {
    plotFilterHistograms(
      dataset = result,
      rawCounts = rawCounts,
      plot.file = file.path(
        parameters$dir$NormalizationImpact,
        paste(sep = "_", parameters$recountID, "filtering_variance_per_gene_hist.pdf")))
  }


  message.with.time("Finished filterDataTable() for Recount experiment ID ", parameters$recountID)
  ################################################################################
  ## Return an object with the filtered counts + updated pheno table selected classes and samples
  ## Transpose the table in order to get it in the suitable format for classifiers:
  ## one row per individual, one column per variable.
  message("\tBefore filtering: ",
          rawCounts$nbGenes, " genes x ",
          rawCounts$nbSamples, " samples",
          " belonging to ", rawCounts$nbClasses, " classes")
  message("\tAfter filtering: ",
          result$nbGenes, " genes x ",
          result$nbSamples, " samples",
          " belonging to ", result$nbClasses, " classes")
  # View(filteredDataTable)
  # dim(filteredDataTable)


  return(result)

}

