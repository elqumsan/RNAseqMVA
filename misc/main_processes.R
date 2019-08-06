library("RNAseqMVA")

###### main steps for the analysis supervised classification methods by RNAseq Data #####
message.with.time("Loading libraries")
source('misc/01a_load_libraries.R')
message.with.time("Loading parameters")
source('misc/01b_load_parameters.R')

if (project.parameters$global$reload) {
  message.with.time("Reloading count data")
  source('misc/02b_reload_counts.R')
} else {
  message.with.time("Loading and normalising raw count data")
  source('misc/02_load_and_normalise_counts.R')
}

## Reload parameters if required (they may have been changed since the
## study cases were built)
if (project.parameters$global$reload.parameters) {
  message.with.time("Reloading parameters")
  source('misc/02b_reload_parameters.R')
}

## Feature selection with first PCs
message.with.time( "Feature selection by first PCs")
source('misc/07_PCA_impact_of_PC_number.R')

## Run analyses with all variables and default parameters
message.with.time( "Impact of normalisation on classifier performances")
source('misc/06_all_variables_vs_all_PCs.R')

# stop("HELLO\tclassifier\t", classifier)


## Test the imapct of kernel on SVM performances
message.with.time( "Impact of kernel on SVM performances")
source('misc/14_svm_impact_of_parameters.R')

## Test the impact of k on KNN performances
message.with.time( "Impact of k on KNN performances")
source('misc/13_knn_impact_of_parameters.R')

message("ALL ANALYSES PERFORMED")


## stop("OLD CODE COMES HEREAFTER")

## if (parameters$compute) {

##   #### Computation of prcomp for the whole raw count #####
##   message.with.time("Computation of prcomp for the raw count data")
##   source("~/RNAseqMVA/misc/03_prcomp_computation.R")
##   message.with.time("Finishing the Computation of prcomp for the raw count data")

##   ##### Computation of DEG (differential expression analysis) for the raw count table #####
##   message.with.time("DEG computation for the raw count table")
##   source("~/RNAseqMVA/misc/04_DEG_computation.R")
##   message.with.time("Finishing from the DEG computation for the raw count table")

##   #### Computation of importance Variables by using random forest ######
##   message.with.time(" computation of the importance variables by RF methods")
##   source("~/RNAseqMVA/misc/05_variables_importance_computation.R")
##   message.with.time("Finishig from computation of the importance variables by RF methods")

##   ##### impact of the number of orderd genes by the edgeR and DESeq2 ####
##   message.with.time("analysis the impact of the number of orderd genes by the edgeR and DESeq2")
##   source("~/RNAseqMVA/misc/08_DEG_impact_of_top_gene_number.R")
##   message.with.time("Finishig analysis the impact of the number of orderd genes by the edgeR and DESeq2")

##   #### impact of the number of importance variables coputed by random forest ####
##   message.with.time(" Analysis the impact of the number of importance variables coputed by random forest")
##   source("~/RNAseqMVA/misc/09_Importance_Variables_impact_of_top_genes.R")
##   message.with.time(" Finishing analysis the impact of the number of importance variables coputed by random forest")



## }
