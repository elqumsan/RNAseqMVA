library("RNAseqMVA")

##### Main steps for the supervised classification of RNA-seq data #####

## Load R libraries
message.with.time("Loading libraries")
source('misc/01a_load_libraries.R')

## Load user-specified parameters for the analysis
message.with.time("Loading parameters")
source('misc/01b_load_parameters.R')

#### Load or reload study cases ####
if (project.parameters$global$reload) {
  message.with.time("Reloading count data")
  source('misc/02b_reload_counts.R')
} else {
  message.with.time("Loading and normalising raw count data")
  source('misc/02_load_and_normalise_counts.R')
}

## Reload parameters if required (they may have been changed since the
## study case was built)
if (project.parameters$global$reload.parameters) {
  message.with.time("Reloading parameters")
  source('misc/01d_reload_parameters.R')
}


## Ensure consistency between iterations and those attached to the study case datasets.
## If they differ, the training / testing sets must be regenerated.
if (studyCase$parameters$iterations < project.parameters$global$iterations) {
  message(studyCase$ID,  "\tResetting the number of iterations from ",
          studyCase$parameters$iterations, " to ", project.parameters$global$iterations)

  ## Set the new number of iterations to the studyCase and to each of its datasets
  studyCase$parameters$iterations <- project.parameters$global$iterations
  for (dataset in studyCase$datasetsForTest) {
    dataset$parameters$iterations <- project.parameters$global$iterations
    ## Regenerate training / testing sets
    buildAttributes(dataset)
  }
}

#### Start parallel computing ####
message.with.time("Initializing parallel computing")
source('misc/01c_init_parallel_computing.R')

#### Feature selection with first PCs ####
message.with.time("Feature selection by first PCs")
source('misc/07_PCA_impact_of_PC_number.R')

#### Test the imapct of kernel on SVM performances ####
message.with.time("Impact of kernel on SVM performances")
source('misc/14_svm_impact_of_parameters.R')


#### Run analyses with all variables and default parameters ####
message.with.time("Impact of normalisation on classifier performances")
source('misc/06_all_variables_vs_all_PCs.R')


#### Test the impact of k on KNN performances ####
message.with.time("Impact of k on KNN performances")
source('misc/13_knn_impact_of_k.R')


message("YEAH! ALL ANALYSES HAVE BEEN PERFORMED")

message("JOB DONE")
