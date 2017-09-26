#' @title  Loading the count table, the Pheno table and vector of classes
#' @author Mustafa ABUELQUMSAN and Jacques van Helden
#' @description Given a ReCount ID, load the count table, the pheno table
#' and the vector of classes, prepares and filters the count table, via
#' wrapper runs by sample.
#' @param recountID  number of experiment
#' @param mergeRuns=T to wrapper the rus by sample
#' @param classColumn = "tissue" ## Column(s) of the pheno table used to define class labels. Passed to filterCountTable()
#' @param minSamplesPerClass=10 min number of samples per class for the filtering.  Passed to filterCountTable()
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
#' @export
loadCounts <- function( recountID,
                        mergeRuns = T,
                        classColumn = "tissue",
                        minSamplesPerClass = 10, ... ) {
  # defineing directories for required libraries
  dir.main <- "~/RNAseqMVA"
  dir.scripts <- file.path(dir.main , "R")
  # loading the required functions
  source(file.path(dir.scripts ,"load_recount_experiment.R" ))
  source(file.path(dir.scripts,"merge_runs.R"))
  source(file.path(dir.scripts,"filterCountTable.R"))

  ################################################
  # loading required libraries and install them if required
  source(file.path(dir.scripts,"required_libraries.R"))
  requiredpackage <- c("caret")
  RequiredCRANPackages(requiredpackage)

  ################################################
  # loading count Data from recount_experiment, Via our wrapper which will Automatically merge the runs by
  # sample in order to obtain sample-wise count rather than run-wise counts.

  countdata <-loadRecountExperiment(recountID=recountID, mergeRuns=mergeRuns, ...)

  ################################################################
  ## Tanspose the count tabe in oredr to use it with classifier methods
  ## which expect a table with one raw per observation (biological samples)
  ## and one column per variable (gene).
  countTable <-as.data.frame(t(countdata$merged$sampleCounts))

  ################################################################
  # Extract pheno table for the sample-wised merged count
  phenoTable <- countdata$merged$samplePheno
  #dim(PhenoTable)

  ################################################
  ## Filter zero-variance and near-zero variance variables from the count table
  filteredData <- filterCountTable(countTable , phenoTable,
                                   classColumn = classColumn, minSamplesPerClass = minSamplesPerClass)

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
  loadedRecount$samples.per.class <- as.data.frame.table(table(classes), row.names=1)
  return(loadedRecount )
}
