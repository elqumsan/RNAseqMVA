#' @title initialise parameters for parallel computing with foreach and doParallel packages
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
