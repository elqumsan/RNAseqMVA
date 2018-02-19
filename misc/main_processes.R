###### main steps for the analysis supervised classification methods by RNAseq Data #####

if(parameters$compute){

  ###### calling the some from miscellaneous file for study the supervised classification methods ####
  ##### loading all parameters and lirraries ####
  message.with.time(" Loading all parameters and libraries")
  source("~/RNAseqMVA/misc/01_parameters_and_libraries.R")
  message.with.time("finishing from Loading all parameters and libraries")

  #### Loading and normalising raw count data ####
  message.with.time("Loading and normalising raw count data")
  source("~/RNAseqMVA/misc/02_load_and_normalise_counts.R")
  message.with.time("Finishing loading and normalising raw count data")

  #### Computation of prcomp for the whole raw count #####
  message.with.time("Computation of prcomp for the raw count data")
  source("~/RNAseqMVA/misc/03_prcomp_computation.R")
  message.with.time("Finishing the Computation of prcomp for the raw count data")

  ##### Computation of DEG (differential expression analysis) for the raw count table #####
  message.with.time("DEG computation for the raw count table")
  source("~/RNAseqMVA/misc/04_DEG_computation.R")
  message.with.time("Finishing from the DEG computation for the raw count table")

  #### Computation of importance Variables by using random forest ######
  message.with.time(" computation of the importance variables by RF methods")
  source("~/RNAseqMVA/misc/05_variables_importance_computation.R")
  message.with.time("Finishig from computation of the importance variables by RF methods")

  #### analysis impact of all variables for the diverse data type onto efficiency of classifiers ####
  message.with.time( "analysis impact of all variables for the diverse data type onto efficiency of classifiers")
  source("~/RNAseqMVA/misc/06_all_variables_vs_all_PCs.R")
  message.with.time( "finishing analysis impact of all variables for the diverse data type onto efficiency of classifiers")

  ###### impact of the number of PCs onto all the classifiers ####
  message.with.time( "analysis the impact of the number of PCs onto the classifiers")
  source("~/RNAseqMVA/misc/07_PCA_impact_of_PC_number.R")
  message.with.time( "finishing analysis the impact of the number of PCs onto the classifiers")

  ##### impact of the number of orderd genes by the edgeR and DESeq2 ####
  message.with.time("analysis the impact of the number of orderd genes by the edgeR and DESeq2")
  source("~/RNAseqMVA/misc/08_DEG_impact_of_top_gene_number.R")
  message.with.time("Finishig analysis the impact of the number of orderd genes by the edgeR and DESeq2")

  #### impact of the number of importance variables coputed by random forest ####
  message.with.time(" Analysis the impact of the number of importance variables coputed by random forest")
  source("~/RNAseqMVA/misc/09_Importance_Variables_impact_of_top_genes.R")
  message.with.time(" Finishing analysis the impact of the number of importance variables coputed by random forest")



}