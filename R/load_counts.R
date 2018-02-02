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
#' @param ...  additional parameters are passed to loadRecountExperiment
#'
#' @return
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

  ################################################
  # loading required libraries and install them if required
  # requiredpackage <- c("caret")
  # RequiredCRANPackages(requiredpackage)

  ################################################
  # loading count Data from recount_experiment, Via our wrapper which will Automatically merge the runs by
  # sample in order to obtain sample-wise count rather than run-wise counts.

  countdata <- loadRecountExperiment(recountID=recountID, mergeRuns=mergeRuns, sampleIdColumn=parameters$sampleIdColumn,
                                     dir.workspace = parameters$dir$workspace, ...)

  ################################################################
  ## Tanspose the count tabe in oredr to use it with classifier methods
  ## which expect a table with one raw per observation (biological samples)
  ## and one column per variable (gene).
  message("Transpose the count table to use it with classifier methods ")
  if(mergeRuns == 1){
  countTable <-as.data.frame(t(countdata$merged$sampleCounts))
  } else {
    countTable <- as.data.frame(t(countdata$countsPerRun))
  }
  #countTable <-as.data.frame(t(countdata$runCounts))


  ################################################################
  # Extract pheno table for the sample-wised merged count
  message("Extracting pheno Table from sample-wised merged count")
  if(mergeRuns == 1){
   phenoTable <- countdata$merged$samplePheno
  } else {
    phenoTable <- countdata$runPhenoTable
  }
  #phenoTable <- countdata$runPheno$characteristics
  countsPerRun <- countdata$countsPerRun
  geo.characteristics <- countdata$geo.characteristics
  #dim(phenoTable)

  ################################################
  ## Filter zero-variance and near-zero variance variables from the count table
  filteredData <- filterCountTable(countTable , phenoTable,
                                   classColumn =parameters$classColumn,
                                   minSamplesPerClass = parameters$minSamplesPerClass)

  # replace unfiltered Data with filted data
  countTable <- filteredData$countTable
  phenoTable <- filteredData$phenoTable

  classes <-make.names(filteredData$classes , unique = F , allow_ = T)
  #classifier requires a factor for classes
  classes <-as.factor(classes)

  message("Count table dimensions: ", nrow(countTable), sep= " x ", ncol(countTable),"\n",
          "Pheno table dimensions: ",nrow(phenoTable), " X ", ncol(phenoTable), "\n",
          "Vector of class labels length: ", length(classes))
  if (length(classes) != nrow(countTable) )
    stop("invaled number of classes labes, (", length(classes) ,"). ",
         "must equal the number of individuals in Count Table (",nrow(countTable),").")

  loadedRecount <- list()
  loadedRecount$countTable<- countTable
  loadedRecount$phenoTable <- phenoTable
  loadedRecount$classes <- classes
  loadedRecount$countsPerRuns <- countsPerRun
  loadedRecount$runPheno <- countdata$runPheno
  loadedRecount$samples.per.class <- as.data.frame.table(table(classes), row.names=1)
  loadedRecount$filteredData <- filteredData
  loadedRecount$geo.characteristics <- geo.characteristics

  message.with.time("Finished Load Count Table process for Recount experiment ID ", parameters$recountID)

  return(loadedRecount )
}
