###### main steps for the analysis supervised classification methods by RNAseq Data #####

if(parameters$compute){

  ###### calling the some from miscellaneous file for study the supervised classification methods ####

  #### analysis impact of all variables for the diverse data type onto efficiency of classifiers ####
  message.with.time( "analysis impact of all variables for the diverse data type onto efficiency of classifiers")
  source("~/RNAseqMVA/misc/06_all_variables_vs_all_PCs.R")

  ###### impact of the number of PCs onto all the classifiers ####
  message.with.time( "impact of the number of PCs onto the classifiers")
  source("~/RNAseqMVA/misc/07_PCA_impact_of_PC_number.R")

  ##### impact of the number of orderd genes by the edgeR and DESeq2 ####
  message.with.time("impact of the number of orderd genes by the edgeR and DESeq2")
  source("~/RNAseqMVA/misc/08_DEG_impact_of_top_gene_number.R")

  #### impact of the number of importance variables coputed by random forest ####
  message.with.time(" impact of the number of importance variables coputed by random forest")
  source("~/RNAseqMVA/misc/09_Importance_Variables_impact_of_top_genes.R")



}
