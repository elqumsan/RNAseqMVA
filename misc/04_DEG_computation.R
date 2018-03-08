
#### Differential analysis with DESeq2 and edgeR to define gene (variable) order ####

if(parameters$compute){
#### Run differential analysis with edgeR to define variable order ####
message.with.time("Running edgeR to define variable ordering")
loaded$DEG.edgeR  <- DEGordering(loaded$originalExperiment$countTable, loaded$originalExperiment$classLabels, method = "edgeR")
message.with.time("Done edgeR DEG analysis")
# sorted.log2.transformed.edgeR <- loaded$countTable[, DEG.DESeq2$geneOrder]

# ######### sptiting the DEG.edgeR dataset for the train set and test set #########
# n <- nrow(DEG.edgeR$orderedCountTable) ## Number of observations (samples)
# train.size <- round(n * parameters$trainingProportion)
#
# ## Random selection of indices for the training set
# trainIndex <- sort(sample(1:n, size=train.size))
# ## Use remaining indices for the testing set
# testIndex <- setdiff(1:n, trainIndex)
#
# DEG.edgeR$trainIndex <- trainIndex
# DEG.edgeR$testIndex <- testIndex

#### Run differential analysis with DESeq2 to define variable order ####
message.with.time("Running DESeq2 to define variable ordering")
loaded$DEG.DESeq2 <- DEGordering(loaded$originalExperiment$countTable, loaded$originalExperiment$classLabels, method = "DESeq2")
message.with.time("Done DESeq2 DEG analysis")
# sorted.log2.transformed.DESeq2 <- loaded$countTable[, DEG.edgeR$geneOrder]

# ######### sptiting the DEG.DESeq2 dataset for the train set and test set #########
# n <- nrow(DEG.DESeq2$orderedCountTable) ## Number of observations (samples)
# train.size <- round(n * parameters$trainingProportion)
#
# ## Random selection of indices for the training set
# trainIndex <- sort(sample(1:n, size=train.size))
# ## Use remaining indices for the testing set
# testIndex <- setdiff(1:n, trainIndex)
#
# DEG.DESeq2$trainIndex <- trainIndex
# DEG.DESeq2$testIndex <- testIndex

#### Order genes at random ####
message.with.time("Generating a fake DEG result by random gene ordering")
DEG.randomized <- list()
## such is randomizing the ordered count table by DEG edgeR; BUT be ware we should keep the neme of this variable as
## DEG.randomized$orderedCountTable to we can correctly use it in the script (DEG_impact_of_top_gene_number.R) where there
## is method to automatically get the name of object under the "$orderedCountTable"; so that you should keep this part
# as it that to avoiding the vey critical error in computation.
loaded$DEG.randomized <- loaded$DEG.edgeR
loaded$DEG.randomized$geneOrder <-sample( loaded$DEG.randomized$geneOrder)
loaded$DEG.randomized$orderedCountTable <- DEG.edgeR$orderedCountTable[ ,sample(DEG.edgeR$geneOrder)]

# ######### sptiting the DEG.edgeR dataset for the train set and test set #########
# n <- nrow(loaded$DEG.randomized$orderedCountTable) ## Number of observations (samples)
# train.size <- round(n * parameters$trainingProportion)
#
# ## Random selection of indices for the training set
# trainIndex <- sort(sample(1:n, size=train.size))
# ## Use remaining indices for the testing set
# testIndex <- setdiff(1:n, trainIndex)
#
# DEG.randomized$trainIndex <- trainIndex
# DEG.randomized$testIndex  <- testIndex

message.with.time("Finishing from the DEG Computation Process")

## define the file to store memory image for the "all DEG Computation" process

image.dir <- file.path (parameters$dir$memoryImages, parameters$recountID)
dir.create(image.dir, showWarnings = FALSE, recursive = TRUE)
image.file <- file.path(image.dir,  paste(sep="","all_DEG_computation_", parameters$recountID,".Rdata"))

if (parameters$save.image) {
  save.image(file = image.file)
}

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

