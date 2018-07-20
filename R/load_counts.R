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
#' @param minSamplesPerClass=parameters$minSamplesPerClass min number of samples per class for the filtering.  Passed to filterDataTable()
#' @param dir.workspace=parameters$dir$workspace Directory in which the data and results will be stored.
#' @param na.rm=TRUE if TRUE, omit all the rows (genes) that contain NA values. True by default because these NA values are problematic for classification methods.
#' @param ...  additional parameters are passed to loadRecountExperiment()
#'
#' @return it is returen the two feature types:
#'  1) for the original counts that means implicitly the raw count table before making ang pre-processing proccess.
#'  2) for the filtered count that means implicitly the raw count after made the pre-processing proccess,
#'
#'   \itemize {
#'      \item filtered: count table after applying the pre-processing proccess.
#'
#'   }
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
#' ## RCurl and XML are Required for recount but not declared in dependencies
#' @import RCurl
#' @import XML
#' @import caret
#' @export
loadCounts <- function(recountID,
                       parameters,
                       ... ) {

  RequiredBioconductorPackages(c("recount", "SummarizedExperiment", "S4Vectors"))

  message.with.time("Starting loadCounts() for Recount ID ", recountID)

  ## Check required parameters
  for (p in c("classColumn", "mergeRuns", "sampleIdColumn", "minSamplesPerClass", "na.rm")) {
    if (is.null(parameters[[p]])) {
      stop("Missing required parameter: '", p,
           "'.\n\tPlease check configuration file. ")
    } else {
      assign(p, parameters[[p]])
    }
  }
  # classColumn = parameters$classColumn
  # mergeRuns = parameters$mergeRuns
  # sampleIdColumn = parameters$sampleIdColumn
  # minSamplesPerClass = parameters$minSamplesPerClass

  ## Check required directory
  if (is.null(parameters$dir$workspace)) {
    stop("Missing required parameter: 'parameters$dir$workspace'.\n\tPlease check configuration file. ")
  } else {
    dir.workspace = parameters$dir$workspace
    dir.create(dir.workspace, recursive = TRUE, showWarnings = FALSE)
  }


  ################################################
  # loading required libraries and install them if required
  # requiredpackage <- c("caret")
  # RequiredCRANPackages(requiredpackage)

  ################################################
  # loading count Data from recount_experiment, Via our wrapper which will Automatically merge the runs by
  # sample in order to obtain sample-wise count rather than run-wise counts.
  experiment <- loadRecountExperiment(recountID = recountID, parameters = parameters, ...)
                                     #  mergeRuns=mergeRuns,
                                     #  sampleIdColumn=sampleIdColumn,
                                     #  dir.workspace = dir.workspace,
                                     #  classColumn = classColumn,
                                     #  na.rm = na.rm,
                                     # ...)


  ## JhV: I suppress this na filtering from here because we want to keep the NA values in
  ## the original count table, and anyway wae do eliminate them in the filtering function.

  # ## If requested, suppress the NA values
  # if (na.rm) {
  #   message("\tSuppressing rows (genes) with NA values")
  #   dataTable <- na.omit(experiment$originalCounts$dataTable)
  #   # dim(dataTable)
  # }

  ## This control is not necessary anymore since this is done when building the DataTableWithClasses object.
  # ##### Check the dimensions of original experiment #####
  # if (experiment$originalCounts$nbSamples != length(experiment$originalCounts$classLabels)){
  #   stop("The number of samples (", experiment$originalCounts$nbSamples, ") differs from the number class labels (", length(experiment$original$classLabels),")")
  # }


  message("\toriginal dataTable contains ",
          nrow(experiment$originalCounts$dataTable), " rows (genes) and ",
          ncol(experiment$originalCounts$dataTable), " columns (samples).")
  #dataTable <-as.data.frame(t(experiment$runCounts))



  ################################################
  #### Filter zero-variance and near-zero variance variables from the count table #####
  experiment$filtered <- filterDataTable(rawCounts = experiment$originalCounts)
  # #,
  #                                         nearZeroVarFilter = parameters$nearZeroVarFilter,
  #                                         minSamplesPerClass = parameters$minSamplesPerClass)
  # class(experiment$filtered)
  # summary(experiment$filtered)
  # dim(dataTable)


  ##### Check the dimensions of filtered experiment #####
  if (length(experiment$filtered$classLabels) != ncol(experiment$filtered$dataTable)) {
    stop("invaled number of filtered class labels, (", length(experiment$filtered$classLabels) ,"). ",
         "must equal the number of individuals in Count Table (",ncol(experiment$filtered$dataTable),").")
  }


  #### Build a list with the results of the loadCounts() function
  experiment$originalCounts$geo.characteristics <- geo_characteristics(experiment$originalCounts$phenoTable)
  experiment$filtered$geo.characteristics <- geo_characteristics(experiment$filtered$phenoTable)





  message.with.time("Finished loadCounts() for Recount experiment ID ", parameters$recountID)
  # experiment$countsPerRun <- NULL
  # experiment$runPhenoTable <- NULL
  # experiment$original <- NULL
  # experiment$filtered <- NULL
  # experiment$samples.per.class <- NULL
  # experiment$merged <- NULL
  # experiment$geo.characteristics <- NULL

  return(experiment)
}
