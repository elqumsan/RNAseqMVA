#' @title initialise parameters for parallel computing with foreach and doParallel packages
#' @author Mustafa AbuELQumsan and Jacques van Helden.
#' @description  by such all required paremater will be initialised for parallel computing by utilizing foreach and doParallel.
#' @param project.parameters that are specified in yaml configuration file in oder to set all default and global variables.
#' @import foreach
#' @import doParallel
#' @return all project.parameters for the related recount ID.
#' \itemize {
#'     \item global: ## Global parameters are those who are the samee across all study cases
#'      \item  dir: ## Directories
#'      \item   main: "~/RNAseqMVA"
#'      \item  workspace: "~/RNAseqMVA_workspace"
#'      \item results: "~/RNAseqMVA_workspace/results"
#'      \item TSV: "~/RNAseqMVA_workspace/data"
#'      \item memoryImages: "~/RNAseqMVA_workspace/memory_images"

#'      \item perm.prefix: "permLabels"
#'      \item jobs: "none" ## Maximal number of jobs to run in parallel. Supported: integers > 0, "auto", or "none"
#  iterations: 50, ## Number of iterations for the classiifers
#'      \item iterations: 10 ## Number of iterations for the classiifers

#'      \item reload: FALSE
#'      \item compute: TRUE ## If FALSE, do not run the heavy computations, just generate the pictures and save tables
#'      \item save.tables: TRUE ## If TRUE, save tab-delimited files with the results
#'      \item save.image: FALSE
#' draw.plot : TRUE
#'      \item classifiers: "svm"
# classifiers: ["svm", "knn","rf"]
# classifiers: ["knn","rf","svm"] # TO ADD LATER: "lda", "qda" + the one suggested by Aitor

#'     \item data.types: ["countsPerRun", "originalCounts", "filtered", "norm", "log2norm", "log2normPCs", "DEG", "V.importance"]
#'      \item data.types.to.test: ["log2normPCs"]
# data.types.to.test: ["filtered", "norm", "log2norm", "log2normPCs"]

# data.types.to.test: ["log2norm"]
#'    \item variables.type: ["all", "top_ordered"]
#'    \item deg.methods: ["DESeq2", "edgeR"]
## Note: for the number of variables I choose a regular spacing around the total number of samples in order to detect overfitting effects
#'    \item nb.variables: [3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 300, 400, 500, 1000, 2000, 5000, 10000] ## Number of variables for variable ordering test
#'    \item trainingProportion: [0.66] # ratio of spliting the data set to training and testing sets
#'    \item identicalTrainTest: TRUE
#'    \item permute: [FALSE, TRUE]
#'    \item verbose: TRUE

#'    \item knn:
#'    \item   k: 10
#'    \item k_values: [3, 5, 7, 10, 15]

#'    \item svm:
#'    \item kernel: "linear"
#'    \item kernel_values: ["linear", "polynomial", "sigmoid"] ## NOT WORKING YET, TO BE TESTED LATER: "radial basis"
#'    \item scale: FALSE
#'    \item type: "C-classification"
# type: ["C-classification", "nu-classification", "one-classification"]
#'    \item gamma: 1
#'    \item cost: 100

## Default values for study case-specific parameters
## (can then be over-written by recountID-specific values as shown below)
#'   \item default:
  #  classColumn: "tissue" ## Characteristics to be used as class label
  #' \item mergeRuns: TRUE ## Whether or not to merge runs per sample
#'   \item sampleIdColumn: "geo_accession"

## Filtering options
#'  \item minSamplesPerClass: 15 ## THERE IS A BUG TO FIX: if we set this parameter to 10 (for example) the class labels must be adapted by suppressing the
#'  \item nearZeroVarFilter: TRUE
#'  \item na.rm: TRUE
#' }
#' @export
initParallelComputing <- function(project.parameters) {
  if (!is.null(project.parameters$global$jobs)) {
    message("\tSetting parallelisation parameters")
    require(foreach)
    require(doParallel)

    ## Define the max number of parallel jobs
    no_cores <- NULL
    jobs <- project.parameters$global$jobs
    if (jobs == "auto") {
      message("\tAutomatic detection of cores for parallel computing. ")
      no_cores <- as.integer(detectCores() - 1) ## Define the number of cores (keep one for the system)
    } else if (jobs == "none") {
      message("\tParallel computing not activated in config file. ")
      project.parameters$global$parallel <- FALSE
    } else {
      jobs <- as.integer (jobs)
      if (is.integer(jobs) && (jobs > 0)) {
        no_cores <- jobs
      } else {
        stop(jobs, ' is not a valid value for the number of jobs. Should be either a strictly positive integer, or "auto".')
      }
    }

    if (is.integer(no_cores) && (no_cores > 0)) {
      message("\tInitializing cluster computing with max ", no_cores, " parallel jobs. ")
      cl <- makeCluster(no_cores) ## Instantiate the cluser
      passed.libs <- clusterEvalQ(cl, library(RNAseqMVA)) ## pass the RNAseqMVA library to the cluster
      clusterExport(cl, "project.parameters") ## pass the global variable "parameters" to the cluster
      registerDoParallel(cl)   ## Register the cluster for parallel computing
      project.parameters$global$cl <- cl
      project.parameters$global$parallel <- TRUE
      project.parameters$global$no_cores <- no_cores
    }

    ## NOTE: WE SHOULD HAVE A TERMINATING SCRIPT, FOR INSTANCE TO STOP THE CLUSTER
    ## Stop the cluster
    #  stopImplicitCluster()

  }
  return(project.parameters)
}
