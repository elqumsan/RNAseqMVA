###### main steps for the analysis supervised classification methods by RNAseq Data #####
message.with.time(" Loading all parameters and libraries")
source('misc/01a_load_libraries.R')
source('misc/01b_load_parameters.R')

message.with.time("Loading and normalising raw count data")
source('misc/02_load_and_normalise_counts.R')

message.with.time( "analysis impact of all variables for the diverse data type onto efficiency of classifiers")
source('misc/06_all_variables_vs_all_PCs.R')

message("ALL ANALYSES PERFORMED")

