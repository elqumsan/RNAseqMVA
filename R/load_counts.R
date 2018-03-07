#' @title  Loading the count table, the Pheno table and vector of class labels
#' @author Mustafa ABUELQUMSAN and Jacques van Helden
#' @description Given a ReCount ID, load the count table, the pheno table
#' and the vector of class labels, prepares and filters the count table, via
#' wrapper runs by sample.
#' @param recountID  number of experiment
#' @param mergeRuns to wrapper the rus by sample
#' @param sampleIdColumn="geo_accession"  name of the column of the pheno table which contains the sample IDs.
#' @param classColumn=parameters$classColumn name of a column of the pheno table which contains the class labels.
#' In some cases classes must be built by concatenating several columns of the pheno table (e.g. "tissue" and "cell.type" for dataset SRP057196),
#' This can be achieved by providing a vector of column names from the pheno table. In this case, class names
## are built by concantenating the values in the specified columns (separated by "_").
#' @param minSamplesPerClass=parameters$minSamplesPerClass min number of samples per class for the filtering.  Passed to filterCountTable()
#' @param dir.workspace=parameters$dir$workspace Directory in which the data and results will be stored.
#' @param na.rm=TRUE if TRUE, omit all the rows (genes) that contain NA values. True by default because these NA values are problematic for classification methods.
#' @param ...  additional parameters are passed to loadRecountExperiment()
#'
#' @return
#'
#' @examples
#'
#' ##############################################
#' ## Load a RNA-seq dataset and merge counts per sample
#' x <- loadCounts( recountID = "SRP048759", mergeRuns=T, classColumn = "tissue")
#'
#' ################################################################
#' ## RNA-seq data combining 2 columns to define sample classes (classLabels vector)
#' ## Reduce the min number of samples per class for this dataset.
#' x <- loadCounts( recountID = "SRP057196", mergeRuns=T,
#'     classColumn = c("tissue", "cell.type"), minSamplesPerClass=5)
#'
#' @import recount
#' @import SummarizedExperiment
#' @import S4Vectors
#' @import caret
#' @export
loadCounts <- function(recountID = parameters$recountID,
                       classColumn = parameters$classColumn,
                       mergeRuns = parameters$mergeRuns,
                       sampleIdColumn = parameters$sampleIdColumn, ## Alternative: use "sample"
                       minSamplesPerClass = parameters$minSamplesPerClass,
                       dir.workspace = parameters$dir$workspace,
                       na.rm = TRUE,
                       ... ) {
  message.with.time("Starting loadCounts() for Recount ID ", recountID)

  ################################################
  # loading required libraries and install them if required
  # requiredpackage <- c("caret")
  # RequiredCRANPackages(requiredpackage)

  ################################################
  # loading count Data from recount_experiment, Via our wrapper which will Automatically merge the runs by
  # sample in order to obtain sample-wise count rather than run-wise counts.
  experiment <- loadRecountExperiment(recountID=recountID,
                                     mergeRuns=mergeRuns,
                                     sampleIdColumn=parameters$sampleIdColumn,
                                     dir.workspace = parameters$dir$workspace,
                                     ...)


  ## If requested, suppress the NA values
  if (na.rm) {
    message("\tSuppressing rows (genes) with NA values")
    countTable <- na.omit(countTable)
    # dim(countTable)
  }


  message("\tcountTable contains ",
          nrow(countTable), " rows (genes) and ",
          ncol(countTable), " columns (samples).")
  #countTable <-as.data.frame(t(experiment$runCounts))



  ################################################
  ## Filter zero-variance and near-zero variance variables from the count table
  experiment$filtered <- filterCountTable(
    countTable = experiment$original$countTable,
    phenoTable = experiment$original$phenoTable,
    classLabels = experiment$original$classLabels,
    #                                   nearZeroVarFilter = FALSE,
    minSamplesPerClass = minSamplesPerClass)
  # dim(countTable)

  # replace unfiltered Data with filted data
  # countTable <- filteredData$countTable
  # phenoTable <- filteredData$phenoTable
  # classLabels <- make.names(filteredData$classLabels , unique = F , allow_ = T)
  # #classifier requires a factor for classLabels
  # classLabels <- as.factor(classLabels) ## RandomForest requires as.factor()

  message("\tUnfiltered count table",
          "\n\t\tdimensions: ", nrow(filteredData$countTable), sep= " x ", ncol(filteredData$countTable),
          "\n\t\tPheno table dimensions: ",nrow(filteredData$phenoTable), " X ", ncol(filteredData$phenoTable),
          "\n\t\tVector of class labels length: ", length(filteredData$classLabels),
          "\n\t\tNumber of distinct classes: ", length(unique(filteredData$classLabels)))
  message("\tFiltered count table",
          "\n\t\tdimensions: ", nrow(filteredData$filteredCountTable), sep= " x ", ncol(filteredData$filteredCountTable),
          "\n\t\tPheno table dimensions: ",nrow(filteredData$filteredPhenoTable), " X ", ncol(filteredData$filteredPhenoTable),
          "\n\t\tVector of class labels length: ", length(filteredData$filteredClasses),
          "\n\tNumber of distinct classes: ", length(unique(filteredData$filteredClasses)))

  if (length(filteredData$filteredClasses) != nrow(filteredData$filteredCountTable)) {
    stop("invaled number of filtered class labels, (", length(filteredData$filteredClasses) ,"). ",
         "must equal the number of individuals in Count Table (",nrow(filteredData$filteredCountTable),").")
  }



  #### Build a list with the results of the loadCounts() function
  #loadedRecount <- list()

  experiment$originalCountTable<- filteredData$countTable
  experiment$originalPhenoTable <- filteredData$phenoTable
  experiment$originalClasses <- filteredData$classLabels

  experiment$filteredCountTable<- filteredData$filteredCountTable
  experiment$filterePhenoTable <- filteredData$filteredPhenoTable
  experiment$filteredClasses <- filteredData$filteredClasses


  experiment$samples.per.class <- as.data.frame.table(table(filteredData$filteredClasses), row.names=1)
  experiment$filteredData <- filteredData
  experiment$geo.characteristics <- geo.characteristics
  experiment$countsPerRuns <- experiment$countsPerRun
  experiment$runPheno <- experiment$runPheno

  message.with.time("Finished loadCounts() for Recount experiment ID ", parameters$recountID)

  return(experiment)
}
