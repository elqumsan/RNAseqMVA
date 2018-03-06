#' @title  Loading the count table, the Pheno table and vector of classes
#' @author Mustafa ABUELQUMSAN and Jacques van Helden
#' @description Given a ReCount ID, load the count table, the pheno table
#' and the vector of classes, prepares and filters the count table, via
#' wrapper runs by sample.
#' @param recountID  number of experiment
#' @param mergeRuns to wrapper the rus by sample
#' @param sampleIdColumn="geo_accession"  name of the column of the pheno table which contains the sample IDs.
#' @param classColumn ## Column(s) of the pheno table used to define class labels. Passed to filterCountTable()
#' @param minSamplesPerClass=parameters$minSamplesPerClass min number of samples per class for the filtering.  Passed to filterCountTable()
#' @param dir.workspace=parameters$dir$workspace Directory in which the data and results will be stored.
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
#' ## RNA-seq data combining 2 columns to define sample classes
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

  ################################################################
  ## Tanspose the count tabe in order to use it with classifier methods
  ## which expect a table with one raw per individual (biological samples)
  ## and one column per feature (gene).
  message("\tTransposing count table to use it with classifier methods. ")
  if (mergeRuns){
    countTable <-as.data.frame(t(na.omit(experiment$merged$sampleCounts)))
  } else {
    countTable <- as.data.frame(t(na.omit(experiment$countsPerRun)))
  }
  message("\tcountTable contains ",
          nrow(countTable), " rows (samples) and ",
          ncol(countTable), " columns (genes).")
  #countTable <-as.data.frame(t(experiment$runCounts))


  ################################################################
  # Extract pheno table for the sample-wised merged count
  message("\tExtracting pheno Table from sample-wised merged count")
  if (mergeRuns){
    phenoTable <- experiment$merged$samplePheno
  } else {
    phenoTable <- experiment$runPhenoTable
  }
  message("\tphenoTable contains ",
          nrow(phenoTable), " rows (samples) and ",
          ncol(phenoTable), " columns (pheno description fields).")
  # names(phenoTable)
  # dim(phenoTable)

  geo.characteristics <- experiment$geo.characteristics
  # dim(geo.characteristics)

  ################################################
  ## Filter zero-variance and near-zero variance variables from the count table
  filteredData <- filterCountTable(countTable = countTable,
                                   phenoTable = phenoTable,
                                   classColumn = classColumn,
#                                   nearZeroVarFilter = FALSE,
                                   minSamplesPerClass = minSamplesPerClass)
  # dim(countTable)

  # replace unfiltered Data with filted data
  # countTable <- filteredData$countTable
  # phenoTable <- filteredData$phenoTable
  # classes <- make.names(filteredData$classes , unique = F , allow_ = T)
  # #classifier requires a factor for classes
  # classes <- as.factor(classes) ## RandomForest requires as.factor()

  message("\tUnfiltered count table",
          "\n\t\tdimensions: ", nrow(filteredData$countTable), sep= " x ", ncol(filteredData$countTable),
          "\n\t\tPheno table dimensions: ",nrow(filteredData$phenoTable), " X ", ncol(filteredData$phenoTable),
          "\n\t\tVector of class labels length: ", length(filteredData$classes),
          "\n\t\tNumber of distinct classes: ", length(unique(filteredData$classes)))
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
  experiment$originalClasses <- filteredData$classes

  experiment$filteredCountTable<- filteredData$filteredCountTable
  experiment$filterePhenoTable <- filteredData$filteredPhenoTable
  experiment$filteredClasses <- filteredData$filteredClasses


  experiment$samples.per.class <- as.data.frame.table(table(filteredData$classes), row.names=1)
  experiment$filteredData <- filteredData
  experiment$geo.characteristics <- geo.characteristics
  experiment$countsPerRuns <- experiment$countsPerRun
  experiment$runPheno <- experiment$runPheno

  message.with.time("Finished Load Count Table process for Recount experiment ID ", parameters$recountID)

  return(experiment)
}
