# ## TEMPOARILY NOT WORKING, HAS TO BE ADAPTED TO TE NEW FUNCTIONS IN R FILES
#
# # self <-filteredDataSet
# # class(self) <- "DataTableWithClasses"
# DataTableWithDEG <- function(self,
#                              DEGmethods = project.parameters$global$ordering.methods,
#                              ...) {
#
#   if (parameters$compute) {
#    # result <-DataTableWithTrainTestSets(self)
#
#     for (method in project.parameters$global$ordering.methods) {
#       DEGdatasets.types[[method]] <- list()
#
#       message.with.time("\tRunning ", method," to define variable ordering","randomized parameter", randomized)
#       DEG.datasets  <- DEGordering(self, method = method, randomized = TRUE)
#       DEGdatasets.types[[method]] <- append(DEGdatasets.types[[method]], DEG.datasets)
#       message.with.time("\t\tDone ",method,"  DEG analysis","randomized parameter", randomized)
#
#     }
#     result$DEGdataSetsType <- DEGdatasets.types
#     result$dataType <- "DEG_order"
#     result$variablesType <- "top_DEG"
#     message.with.time("Finishing from Processing for instantiate DEG object")
#     class(result) <- unique(c("DataTableWithDEG","DataTableWithTrainTestSets"))
#     return(result)
#
#
#   } # end of computation
# } # end of the constructor function
#
# #  if (method == "edgeR"){
# #### Run differential analysis with edgeR to define variable order ####
# message.with.time("\tRunning edgeR to define variable ordering")
# result$DEG.edgeR  <- DEGordering(self$dataTable, self$classLabels, method = "edgeR")
#
# # result$DEG.edgeR$method <-"edgeR"
# message.with.time("\t\tDone edgeR DEG analysis")
# # sorted.log2.transformed.edgeR <- studyCases$dataTable[, DEG.DESeq2$geneOrder]
#
# #### Order genes at random ####
# message.with.time("\tGenerating a fake DEG.edgeR result by random gene ordering")
# DEG.edgeR.randomized <- list()
# ## such is randomizing the ordered count table by DEG edgeR; BUT be ware we should keep the neme of this variable as
# ## DEG.randomized$orderedDataTable to we can correctly use it in the script (DEG_impact_of_top_gene_number.R) where there
# ## is method to automatically get the name of object under the "$orderedDataTable"; so that you should keep this part
# # as it that to avoiding the vey critical error in computation.
#
# # result$DEG.edgeR.randomized <- result$DEG.edgeR
# # result$DEG.edgeR.randomized$geneOrder <-sample( result$DEG.edgeR.randomized$geneOrder)
# # result$DEG.edgeR.randomized$orderedDataTable <- result$DEG.edgeR$orderedDataTable[ sample(result$DEG.edgeR$geneOrder),]
# # result$DEG.edgeR.randomized$method <- "randomaised-edgeR"
# # } else {
#
#
# #   #### Run differential analysis with DESeq2 to define variable order ####
# #   message.with.time("\tRunning DESeq2 to define variable ordering")
# #   result$DEG.DESeq2 <- DEGordering(self$dataTable, self$classLabels, method = "DESeq2")
# # # result$DEG.DESeq2$method <- "DESeq2"
# #   message.with.time("\t\tDone DESeq2 DEG analysis")
# #   # sorted.log2.transformed.DESeq2 <- studyCases$dataTable[, DEG.edgeR$geneOrder]
# #
# #   message.with.time("\tGenerating a fake DEG.DESeq2 result by random gene ordering")
# #   DEG.DESeq2.randomized <- list()
# ## such is randomizing the ordered count table by DEG edgeR; BUT be ware we should keep the neme of this variable as
# ## DEG.randomized$orderedDataTable to we can correctly use it in the script (DEG_impact_of_top_gene_number.R) where there
# ## is method to automatically get the name of object under the "$orderedDataTable"; so that you should keep this part
# # as it that to avoiding the vey critical error in computation.
#
#
# # result$DEG.DESeq2.randomized <- result$DEG.DESeq2
# # result$DEG.DESeq2.randomized$geneOrder <-sample( result$DEG.DESeq2.randomized$geneOrder)
# # result$DEG.DESeq2.randomized$orderedDataTable <- result$DEG.DESeq2$orderedDataTable[sample(result$DEG.DESeq2.randomized$geneOrder) ,]
# # result$DEG.DESeq2.randomized$method < "randomized-DESeq2"
# # }
# # result$DEG.datasets<- list(result$DEG.edgeR, result$DEG.edgeR.randomized$orderedDataTable, result$DEG.DESeq2, result$DEG.DESeq2.randomized$orderedDataTable)
# # names(result$DEG.datasets) <- c("DEG.edgeR", "DEG.edgeR.randomized","DEG.DESeq2","DEG.DESeq2.randomized")
#
#
#
# ## define the file to store memory image for the "all DEG Computation" process
#
# image.dir <- file.path (parameters$dir$memoryImages, parameters$recountID)
# dir.create(image.dir, showWarnings = FALSE, recursive = TRUE)
# image.file <- file.path(image.dir,  paste(sep="","all_DEG_computation_", parameters$recountID,".Rdata"))
#
# if (parameters$save.image) {
#   save.image(file = image.file)
#
#   ##### if compution not required, you can load the image file without any computations ####
# } else {
#   # reload previous results if exist
#   if (file.exists(image.file)) {
#     message ("Reloading memory image ", image.file)
#     load(image.file)
#   } else {
#     stop("Cannot reload memory image file ", image.file)
#   }
# }
#
