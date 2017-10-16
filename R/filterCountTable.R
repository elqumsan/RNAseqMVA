#' @title Apply various filters on the observations and variables of
#' a count table in order to prepare it for classification.
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @description Apply various filters on the observations (e.g. biological samples)
#' and variables (e.g. genes) of
#' a count table in order to prepare it for classification.
#' @param countTable a data.frame with one row per observation and one column per variable
#' @param phenoTable a data.frame with one row per observation and one column per attribute
#' @param classColumn="tissue" name of a column of the pheno table which contains the class labels (default: "tissue").
#' In som cases classes must be built by concatenating several columns of the pheno table (e.g. "tissue" and "cell.type" for dataset SRP057196),
#' This can be achieved by providing a vector of column names from the pheno table. In this case, class names are built by
## are built by concantenating the values in the specified columns (separated by "_").
#' @param minSamplesPerClass=10 minimum nuimber of samples per class to keep
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
#' @export
filterCountTable <- function(countTable, phenoTable,
                             classColumn="tissue",
                             minSamplesPerClass = 10, ...) {

  ## Check if there are NA values, and discard all genes having at least one NA value
  if (sum(is.na(countTable)) > 0) {
    naGenes <- colnames(countTable)[apply(is.na(countTable), 2, sum) > 0]
    message("Filtering out ", length(naGenes)," genes with NA values")
    if (naGenes > 0) {
      keptGenes <- setdiff(colnames(countTable), naGenes)
      countTable <- countTable[, keptGenes]
      dim(countTable)
    }
  }

  ## Count the number of genes with zero variance
  countvar <- apply(countTable, 2 , var, na.rm=TRUE) # Compute variance per gene (row)
  # sum(countvar == 0) # count genes with zero variance
  zeroVarGenes <- names(countvar)[countvar == 0] # identify genes having zero variance
  length(zeroVarGenes) ## Number of genes to discard
  keptGenes <- names(countvar)[countvar > 0]
  #length(keptGenes) ## Number of genes to keep
  countTable <- countTable[, keptGenes]
  #dim(countTable)
  message("Filtering out ", length(zeroVarGenes)," genes with zero variance; keeping ", length(keptGenes), " genes")

  ## Use caret::nearZeroVar() to discard supposedly bad predictor genes.
  ## This includes zero variance genes (already filtered above) but also less trivial cases, see nearZeroVar() doc for details.
  nearZeroVarGenes <- colnames(countTable)[nearZeroVar(countTable)]
  #length(nearZeroVarGenes)
  keptGenes <- setdiff(colnames(countTable), nearZeroVarGenes)
  message("Filtering out ", length(nearZeroVarGenes)," poor predictor genes with caret::nearZeroVar(); kept genes: ", length(keptGenes))
  countTable <- countTable[, keptGenes]
  # dim(countTable)


  sample.nb <- nrow(countTable)
  gene.nb <- ncol(countTable)
  # dim(countTable)
  # dim(phenoTable)
  # View( phenoTable)


  # Specify sample classes by extracting information about specified class columns
  if (length(classColumn) == 1) {
    classes <-  phenoTable[, classColumn]
  } else if (length(classColumn) > 1) {
    ## Combine several columns to establis the classes
    classes <- apply(phenoTable[, classColumn], 1, paste, collapse="_")
  }
  # table(classes)


  ################################################################
  ## Check if there are NA values in the sample classes
  discardedSamples <- is.na(classes)
  message("Discarding ", sum(discardedSamples), " samples with undefined class in ", classColumn, " column of pheno table")
  if (sum(discardedSamples) > 0) {
    nonaSamples <- !is.na(tissue)
    countTable <- countTable[nonaSamples, ]
    phenoTable<- phenoTable[nonaSamples, ]
    classes <- classes[nonaSamples]
    sample.nb <- nrow(countTable)
  }

  ################################################################
  ## Select classess for which we dispose of at least 10 samples
  # length(classes)
  message("Selecting classes with at least ", minSamplesPerClass, " samples")
  samplesPerClass <- table(classes)
  selectedClasses <- names(samplesPerClass)[samplesPerClass >= minSamplesPerClass]
  selectedSamples <- classes %in% selectedClasses
  message("Keeping ",
          sum(selectedSamples), " samples from ",
          length(selectedClasses), " classes (",
          paste(selectedClasses, collapse=", "), ")")

  # update count table, pheno table and classes vector to keep only the samples belonging to selected classess
  countTable <- countTable[selectedSamples, ]
  phenoTable<- phenoTable[selectedSamples, ]
  classes <- classes[selectedSamples]
  sample.nb <- nrow(countTable)

  ################################################################################
  ## Build a countTable with the selected samples.
  ## Transpose the table in order to get it in the suitable format for classifiers:
  ## one row per individual, one column per variable.
  message("After filtering, count table contains ",
          nrow(countTable), " samples (rows) and ",
          ncol(countTable), " genes (columns) ",
          "belonging to ", length(selectedClasses), " classes")
  # View(countTable)
  # dim(countTable)

  result <- list()
  result$countTable <- countTable
  result$phenoTable <- phenoTable
  result$classes <- classes
  result$selectedClasses <- selectedClasses
  return(result)

}
