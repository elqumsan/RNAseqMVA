#################### Load counts and pheno ####################
## Load a count Table from recount-experiment, merge counts per sample
## and apply some pre-filtering (remove zero-variance and near-zero-variance genes).
if (parameters$compute) {
  message.with.time("Loading count table from recount", "; recountID = ", parameters$recountID)
  loaded <- loadCounts(recountID = parameters$recountID,
                       mergeRuns = TRUE,
                       classColumn = parameters$classColumn,
                       minSamplesPerClass = parameters$minSamplesPerClass,
                       na.rm = parameters$na.rm)



  ## Class colors may be defined in the yaml parameters
  if (!is.null(parameters$classColor)) {
    #  parameters$classColor[["astrocytes"]]
    ## Convert the yaml-imported list into a named vector
    loaded$filtered$classColors <- unlist(parameters$classColor) # convert list to a vector
    # names(classColors) # check vector names
    # classColors["astrocytes"]

  } else {
    classNames <- unique(classLabels)
    classColors <- 1:length(classNames)
    names(classColors) <- classNames
  }


  # dim(loaded$originalCountTable)
  # dim(loaded$filteredCountTable)

  ############################################################
  ## Build an object (formally a simple list) for the raw counts table and associated info.
  ## We will then treat similar objects for different types of pre-processed data:
  ## - log2-transformed
  ## - normalised
  ## - CPA-transformed
  ## - ....

  ## JvH 2018-03-06: this rawCounts list is redundant with loaded$filteredCountTable
  # rawCounts <- list()
  # rawCounts$Counts <- loaded$filteredCountTable ## Note: one row per sample, one column per gene
  # rawCounts$sample.nb <- nrow(rawCounts$Counts)
  # rawCounts$feature.nb <- ncol(rawCounts$Counts)
  # dim(rawCounts$Counts)


  ## Assign a specific color to each sammple according to its class
  countTable <- loaded$filtered$countTable
  pheno <- loaded$filtered$phenoTable
  classes <- loaded$filtered$classLabels
  geo.characteristics <- loaded$filtered$geo.characteristics
  # table(classes)
  # length(classes)
  # length(unique(classes))
  distinct.classes <- as.vector(unique(classes))




} else {
  message.with.time("Skipping data loading")
}



## Number of samples per class
# print(loaded$filtered$nbSamples)

## I start by assigning one systematic color(number) to each class,
## To make sure that each class has a color even with different datasets analsed in the efuture.
classColors <- 1:length(distinct.classes)
names(classColors) <- distinct.classes

## JvH: Mustafa, these colors are specific for one dataset.
##
## Message from JvH, 2018-03-06:
## These should now be placed in the yaml file.
if (parameters$recountID == "SRP048759") {
  classColors["Heparinised.blood"] <- "#BB0000"
  classColors["Bone.marrow"] <- "#4488FF"
} else if (parameters$recountID == "SRP042620") {
  ## TO DO : define colors for the multi-group breast cancer dataset
  classColors["Breast.Cancer.Cell.Line"] <- "red"
  classColors["ER..Breast.Cancer.Primary.Tumor"] <- "darkblue"
  classColors["Triple.Negative.Breast.Cancer.Primary.Tumor"] <- "purple"
  classColors["Uninvolved.Breast.Tissue.Adjacent.to.ER..Primary.Tumor"] <- "green"
  classColors["Uninvolved.Breast.Tissue.Adjacent.to.TNBC.Primary.Tumor"] <- "black"
}
# print(classColors)
## Assign colors per sample according to their class
sampleColors <- classColors[loaded$filtered$classLabels]
names(sampleColors) <- rownames(loaded$filtered$phenoTable)
# print(sampleColors)

##### Normalize the counts without log2 transformation (first test) #####
##
## Note: this method takes a table with one column per sample and one
## row per gene, we thus have to transpose the raw count table.
if (parameters$compute) {
  # dim(loaded$loaded$norm$counts)
  message.with.time("Normalizing counts based on 75th percentile")
  loaded$norm <- NormalizeCounts(
    objectFiltered = loaded$filtered,
    # phenoTable = loaded$filteredExperiment$phenoTable,
    # classLabels = loaded$filteredExperiment$classLabels,
    method = "quantile", quantile=0.75, log2 = FALSE)
  # dim(loaded$normCounts)
  # loaded$norm$nb.samples
  # loaded$norm$nb.genes


#  hist(unlist(loaded$loaded$norm$counts), main="Normalised count distribution", breaks=1000)

} else {
  message.with.time("Skipping normalisation")
}



##### Normalize counts with log2 transformation (second test) #####
##
## Note: this method takes a table with one column per sample and one
## row per gene, we thus have to transpose the raw count table.
if (parameters$compute) {
  message.with.time("Normalizing counts based on 75th percentile + log2 transformation")
  loaded$log2norm <- NormalizeCounts(
    objectFiltered = loaded$filtered,
    # counts =loaded$filteredExperiment$countTable,
    # phenoTable = loaded$filteredExperiment$phenoTable,
    # classLabels = loaded$filteredExperiment$classLabels,
    method = "quantile", quantile=0.75,
    log2 = TRUE, epsilon=0.1)
  dim(loaded$log2norm$counts)




  # plot.file <- file.path(dir.NormImpact, "log2normCount_hist.pdf")
  # message("\tlog2(norm counts) histogram\t", plot.file)
  # pdf(plot.file, width=7, height=5)
  # hist(unlist(loaded$log2norm$counts), breaks=100,
  #      col="grey",
  #      main=paste("log2(norm counts) distrib;", recountID),
  #      las=1,
  #      xlab="log2(norm counts)",
  #      ylab="Frequency")
  #
  # silence <- dev.off()

  if (ncol(loaded$log2norm$counts) != length(loaded$log2norm$classLabels)){
    stop(" the Number of samples in log2norm counts should be the same length of classes")
  }



} else {
  message.with.time("Skipping normalisation with log2 transformation")
}



## Indicate that this script has finished running
message.with.time("finished executing 02_load_and_normalise_counts.R")

##### plotting some figures to explore the nuture of recount data set #####
# message.with.time(" plotting some figures to explore distribution for the recount data set ",parameters$recountID)
# source("misc/11_impact_of_normalization_and_Log2.R")

##### Exhibiting the geo charactiristics for the current project #####
message.with.time("Exhibit the geo charactiristics for such experiment: ", parameters$recountID, " in order to know the class lable
                  for such experiment")
head( geo.characteristics)
geo.characteristics.file <- file.path("~/RNAseqMVA_workspace", "data", parameters$recountID, "geo.characteristics.tsv")
write.table( geo.characteristics, file = geo.characteristics.file, quote = FALSE,
             row.names = FALSE, sep = "\t" )
