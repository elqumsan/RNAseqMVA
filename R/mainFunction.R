####################################################################
#' @title main script to upload all required function that are necessary in RNAseqMVA package
#' @author Mustafa ABUELQUMSAN and Jacques van-Helden
#' @description we would like to ensure all function are uploaded in proper way to be aware all analysis steps will
#' # work explicitly way.
#' @export
#'
#'

if (exists( "DEGordering" ,  mode="function")){

  ## This is no required anymore since functions are loaded via the R package
  # source("~/RNAseqMVA/R/deg_ordering.R")
  # source("~/RNAseqMVA/R/ErrorRateBoxPlot.R")
  # source("~/RNAseqMVA/R/load_counts.R")
  # source("~/RNAseqMVA/R/merge_runs.R")
  # source("~/RNAseqMVA/R/misclassification_estimate.R")
  # source("~/RNAseqMVA/R/normalize_count_table.R")
  # source("~/RNAseqMVA/R/normalize_count_table.R")
  # source("~/RNAseqMVA/R/IterateTrainingTesting.R")
  # source("~/RNAseqMVA/R/required_libraries.R")
  # source("~/RNAseqMVA/R/required_libraries.R")

#
# } else if(exists("ErrorRateBoxPlot" , mode = "function")){
#   source("~/RNAseqMVA/R/ErrorRateBoxPlot.R")
# } else if (exists("loadCounts" , mode = "function")){
#   source("~/RNAseqMVA/R/load_counts.R")
# }else if( exists("loadRecountExperiment", mode = "function")){
#     source("~/RNAseqMVA/R/load_recount_experiment.R")
# }else if ( exists("MergeRuns" , mode = "function")){
#     source("~/RNAseqMVA/R/merge_runs.R")
# }else if ( exists("MisclassificationEstimate", mode = "function")){
#     source("~/RNAseqMVA/R/misclassification_estimate.R")
# }else if( exists("NormalizeSamples" ,mode = "function")){
#     source("~/RNAseqMVA/R/normalize_count_table.R")
# }else if( exists("IterateTrainingTesting" , mode = "function")){
#     source("~/RNAseqMVA/R/IterateTrainingTesting.R")
# }else if(exists("RequiredCRANPackages", mode = "function")){
#     source("~/RNAseqMVA/R/required_libraries.R")
# } else if (exists("RequiredBioconductorPackages" , mode = "function")){
#     source("~/RNAseqMVA/R/required_libraries.R")
}else {
    message("such function are not provided in the RNAseqMVA package,!!!")
}
