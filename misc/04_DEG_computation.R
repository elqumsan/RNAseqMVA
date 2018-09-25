## This script takes as input a pre-loaded study cases
## and runs differential expression on the filtered counts
##
## Authors: Mustafa AbuElQumsan and Jacques van Helden
##
## Prerequisite: study cases can be loaded either by running
## the scripts 02_load_and_normalise_counts.R, or by reloading
## a memory image previously stored by this script.
##
## Example
##   load("../RNAseqMVA_workspace/memory_images/loaded_studyCases_SRP042620_[DATE].Rdata")
##    recountID <- "SRP042620"

## TO DO
## - Mustafa: update the flow chart of the scripts,
##   and show the connection between this one and 02

## The flowchart shows the chaining between operations.
## The class diagram shows the oganisation of the objects: methods, attribtues and links between classes.

## Get the selected study case
studyCase <- studyCases[[recountID]]
# class(studyCase) ## Check the class (shoud be StudyCase)
# attributes(studyCase) ## Check the attributes

## Select the filtered dataset for this study case
filteredDataSet <- studyCase["datasetsForTest"]$datasetsForTest[["filtered"]]
# class(filteredDataSet)
# attributes(filteredDataSet)

## Check the DEG methods specified in the YAML configuration file
deg.methods <- project.parameters$global$deg_analysis$methods
if (is.null(deg.methods)) {
  stop("At least one DEG method should be specified in the YAML configuration file")
}
message("\tDEG method(s)\t", paste(collapse=", ", deg.methods))

## Check that alpha was defined, if not use a standard default value
alpha <- project.parameters$global$deg_analysis$alpha
if (is.null(alpha)) {
  alpha <- 0.05
  message("The alpha parameter was not defined in the YAML config file, using default. ")
}
message("\talpha\t", alpha)


deg.results <- list()

#### edgeR-sorted variables ####
for (method in deg.methods) {
  message.with.time("Running differential analysis with ", method)
  deg.results[[method]] <- DEGordering(dataTableWithClasses = filteredDataSet, method = method)
  # class(deg.results[[method]])
  # names(deg.results[[method]])
  # View(deg.results[[method]])


  hist(deg.results[[method]]$DEGtable$padj, breaks=seq(0, 1, 0.05))
  hist(deg.results[[method]]$DEGtable$log2FC, breaks = 100)

  ## Select the genes declared positive (null hypothesis rejected)
  ## Select genes that pass the alpha threshold on adjuste p-value
  alpha.selection <- deg.results[[method]]$DEGtable$padj < alpha
  sum(alpha.selection)

  ## Select genes passing the threshold on log2 fold change
  ## TO DO

  ## Gene clustering

  ## Draw a heatmap of clustered genes

  ## Sample clustering (with only the genes declared positive)
}




#### Differential analysis with DESeq2 and edgeR to define gene (variable) order ####
DataTableWithDEG <- function(self,
                              DEGmethods = project.parameters$global$ordering.methods,
                              ...){

  if(parameters$compute){
    result <-DataTableWithTrainTestSets(self)

    for(method in DEGmethods){
      DEGdatasets.types[[method]] <- list()

      message.with.time("\tRunning ", method," to define variable ordering","randomized parameter", randomized)
      DEG.datasets  <- DEGordering(self$dataTable, self$classLabels, method = method, randomized = TRUE )
      DEGdatasets.types[[method]] <- append(DEGdatasets.types[[method]], DEG.datasets  )
      message.with.time("\t\tDone ",method,"  DEG analysis","randomized parameter", randomized)

    }
    result$DEGdataSetsType <- DEGdatasets.types
    result$dataType <- "DEG_order"
    result$variablesType <- "top_DEG"
    message.with.time("Finishing from Processing for instantiate DEG object")
    class(result) <- unique(c("DataTableWithDEG","DataTableWithTrainTestSets"))
    return(result)


  } # end of computation
} # end of the constructor function

      #  if (method == "edgeR"){
      #### Run differential analysis with edgeR to define variable order ####
      message.with.time("\tRunning edgeR to define variable ordering")
      result$DEG.edgeR  <- DEGordering(self$dataTable, self$classLabels, method = "edgeR")

      # result$DEG.edgeR$method <-"edgeR"
      message.with.time("\t\tDone edgeR DEG analysis")
      # sorted.log2.transformed.edgeR <- studyCases$dataTable[, DEG.DESeq2$geneOrder]

      #### Order genes at random ####
      message.with.time("\tGenerating a fake DEG.edgeR result by random gene ordering")
      DEG.edgeR.randomized <- list()
      ## such is randomizing the ordered count table by DEG edgeR; BUT be ware we should keep the neme of this variable as
      ## DEG.randomized$orderedDataTable to we can correctly use it in the script (DEG_impact_of_top_gene_number.R) where there
      ## is method to automatically get the name of object under the "$orderedDataTable"; so that you should keep this part
      # as it that to avoiding the vey critical error in computation.

      # result$DEG.edgeR.randomized <- result$DEG.edgeR
      # result$DEG.edgeR.randomized$geneOrder <-sample( result$DEG.edgeR.randomized$geneOrder)
      # result$DEG.edgeR.randomized$orderedDataTable <- result$DEG.edgeR$orderedDataTable[ sample(result$DEG.edgeR$geneOrder),]
      # result$DEG.edgeR.randomized$method <- "randomaised-edgeR"
    # } else {


    #   #### Run differential analysis with DESeq2 to define variable order ####
    #   message.with.time("\tRunning DESeq2 to define variable ordering")
    #   result$DEG.DESeq2 <- DEGordering(self$dataTable, self$classLabels, method = "DESeq2")
    # # result$DEG.DESeq2$method <- "DESeq2"
    #   message.with.time("\t\tDone DESeq2 DEG analysis")
    #   # sorted.log2.transformed.DESeq2 <- studyCases$dataTable[, DEG.edgeR$geneOrder]
    #
    #   message.with.time("\tGenerating a fake DEG.DESeq2 result by random gene ordering")
    #   DEG.DESeq2.randomized <- list()
      ## such is randomizing the ordered count table by DEG edgeR; BUT be ware we should keep the neme of this variable as
      ## DEG.randomized$orderedDataTable to we can correctly use it in the script (DEG_impact_of_top_gene_number.R) where there
      ## is method to automatically get the name of object under the "$orderedDataTable"; so that you should keep this part
      # as it that to avoiding the vey critical error in computation.


      # result$DEG.DESeq2.randomized <- result$DEG.DESeq2
      # result$DEG.DESeq2.randomized$geneOrder <-sample( result$DEG.DESeq2.randomized$geneOrder)
      # result$DEG.DESeq2.randomized$orderedDataTable <- result$DEG.DESeq2$orderedDataTable[sample(result$DEG.DESeq2.randomized$geneOrder) ,]
      # result$DEG.DESeq2.randomized$method < "randomized-DESeq2"
   # }
    # result$DEG.datasets<- list(result$DEG.edgeR, result$DEG.edgeR.randomized$orderedDataTable, result$DEG.DESeq2, result$DEG.DESeq2.randomized$orderedDataTable)
    # names(result$DEG.datasets) <- c("DEG.edgeR", "DEG.edgeR.randomized","DEG.DESeq2","DEG.DESeq2.randomized")



## define the file to store memory image for the "all DEG Computation" process

image.dir <- file.path (parameters$dir$memoryImages, parameters$recountID)
dir.create(image.dir, showWarnings = FALSE, recursive = TRUE)
image.file <- file.path(image.dir,  paste(sep="","all_DEG_computation_", parameters$recountID,".Rdata"))

if (parameters$save.image) {
  save.image(file = image.file)

  ##### if compution not required, you can load the image file without any computations ####
} else {
  # reload previous results if exist
  if (file.exists(image.file)) {
    message ("Reloading memory image ", image.file)
    load(image.file)
  } else {
    stop("Cannot reload memory image file ", image.file)
  }
}

