


#### Initialise cluster processing ####
if (!is.null(project.parameters$global$jobs)) {
  ## NOTE: WE SHOULD HAVE A TERMINATING SCRIPT, FOR INSTANCE TO STOP THE CLUSTER
  ## Stop the cluster
  #  stopImplicitCluster()
  project.parameters <- initParallelComputing(project.parameters)
  cl <- project.parameters$global$cl
}

