#' @title ordaring the varaibles by using the P-Value for all variables
#' @author Mustafa ABUELQUMSAN and Jacques van Helden
#' @description  in such constructor we build count table with ordaring the variables based on the most significance.
#' @param self it is filtered count table object
#'
#'
#' @export
#### Differential analysis with DESeq2 and edgeR to define gene (variable) order ####
countTableWithDEG <- function(self,
                              method = "edgeR"){

  if(parameters$compute){
    result <-list()

    if (method == "edgeR"){
      #### Run differential analysis with edgeR to define variable order ####
      message.with.time("Running edgeR to define variable ordering")
      result$DEG.edgeR  <- DEGordering(self$countTable, self$classLabels, method = method)
      message.with.time("\t\tDone edgeR DEG analysis")
      # sorted.log2.transformed.edgeR <- loaded$countTable[, DEG.DESeq2$geneOrder]

      #### Order genes at random ####
      message.with.time("Generating a fake DEG result by random gene ordering")
      DEG.edgeR.randomized <- list()
      ## such is randomizing the ordered count table by DEG edgeR; BUT be ware we should keep the neme of this variable as
      ## DEG.randomized$orderedCountTable to we can correctly use it in the script (DEG_impact_of_top_gene_number.R) where there
      ## is method to automatically get the name of object under the "$orderedCountTable"; so that you should keep this part
      # as it that to avoiding the vey critical error in computation.
      result$DEG.edgeR.randomized <- result$DEG.edgeR
      result$DEG.edgeR.randomized$geneOrder <-sample( result$DEG.edgeR.randomized$geneOrder)
      result$DEG.edgeR.randomized$orderedCountTable <- result$DEG.edgeR$orderedCountTable[ sample(result$DEG.edgeR$geneOrder),]

    } else {


      #### Run differential analysis with DESeq2 to define variable order ####
      message.with.time("Running DESeq2 to define variable ordering")
      result$DEG.DESeq2 <- DEGordering(self$countTable, self$classLabels, method = method)
      message.with.time("\t\tDone DESeq2 DEG analysis")
      # sorted.log2.transformed.DESeq2 <- loaded$countTable[, DEG.edgeR$geneOrder]

      message.with.time("Generating a fake DEG result by random gene ordering")
      DEG.DESeq2.randomized <- list()
      ## such is randomizing the ordered count table by DEG edgeR; BUT be ware we should keep the neme of this variable as
      ## DEG.randomized$orderedCountTable to we can correctly use it in the script (DEG_impact_of_top_gene_number.R) where there
      ## is method to automatically get the name of object under the "$orderedCountTable"; so that you should keep this part
      # as it that to avoiding the vey critical error in computation.
      result$DEG.DESeq2.randomized <- result$DEG.DESeq2
      result$DEG.DESeq2.randomized$geneOrder <-sample( result$DEG.DESeq2.randomized$geneOrder)
      result$DEG.DESeq2.randomized$orderedCountTable <- result$DEG.DESeq2$orderedCountTable[sample(result$DEG.DESeq2.randomized$geneOrder) ,]

    }

    message.with.time("Finishing from the DEG Computation Process")

  }
} # end of the constructor function

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

