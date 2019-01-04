#' @title  Loading the count table, the Pheno table and vector of class labels
#' @author Mustafa ABUELQUMSAN and Jacques van Helden
#' @description Given a ReCount ID, load the count table, the pheno table
#' and the vector of class labels, prepares and filters the count table, via
#' wrapper runs by sample.
#' @param recountID  number of experiment
#' @param parameters list of parameters (supposedly loaded from the yaml configuration file)
#' @param ...  additional parameters are passed to loadRecountExperiment()
#'
#' @return it is returen the two feature types:
#'  1) for the original counts that means implicitly the raw count table before making ang pre-processing proccess.
#'  2) for the filtered count that means implicitly the raw count after made the pre-processing proccess,
#'
#'   \itemize{
#'      \item filtered: count table after applying the pre-processing proccess.
#'
#'   }
#' @examples
#'
#' ##############################################
#' ## Load an RNA-seq dataset and merge counts per sample.
#' This assumes that the required parameters have been specified
#' in the YAML configuration file.
#' x <- loadCounts( recountID = "SRP048759", parameters = parameters)
#'
#'
#' ## RCurl and XML are Required for recount but not declared in dependencies
#' @import RCurl
#' @import XML
#' @import caret
#' @export
loadCounts <- function(recountID,
                       parameters,
                       ... ) {

  LoadRequiredBioconductorPackages(c("recount", "SummarizedExperiment", "S4Vectors"))

  message.with.time("Starting loadCounts() for Recount ID ", recountID)

  ## Check required parameters
  for (p in c("classColumn", "mergeRuns", "sampleIdColumn", "feature")) {
    if (is.null(parameters[[p]])) {
      stop("Missing required parameter: '", p,
           "'.\n\tPlease check configuration file. ")
    } else {
      assign(p, parameters[[p]])
    }
  }

  ## Check required directory
  if (is.null(parameters$dir$workspace)) {
    stop("Missing required parameter: 'parameters$dir$workspace'.\n\tPlease check configuration file. ")
  } else {
    dir.workspace = parameters$dir$workspace
    dir.create(dir.workspace, recursive = TRUE, showWarnings = FALSE)
  }


  ################################################
  # loading count Data from recount_experiment, Via our wrapper which will Automatically merge the runs by
  # sample in order to obtain sample-wise count rather than run-wise counts.
  experiment <- loadRecountExperiment(recountID = recountID, parameters = parameters, ...)

  message("\t", "Original dataTable contains ",
          nrow(experiment$originalCounts$dataTable), " rows (", parameters$feature,") and ",
          ncol(experiment$originalCounts$dataTable), " columns (samples).")

  ################################################
  #### Filter zero-variance and near-zero variance variables from the count table #####
  experiment$filtered <- filterDataTable(rawCounts = experiment$originalCounts) #, parameters)


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
