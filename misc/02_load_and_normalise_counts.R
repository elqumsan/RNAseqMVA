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


  ## Select training and testing sets on the filtered table with raw counts
  ## These wil then be passed to all the derived count tables (normalised, DGE, ...)
  loaded$filtered <- countTableWithTrainTestSets(loaded$filtered)

  #### Export the count tables with their associated information (pheno table, class labels) in tab-separated value (.tsv) files ###
  exportTables(loaded$countsPerRun,
               export.dir = paste(parameters$dir$TSV, parameters$recountID, sep = "/"),
               file.prefix = "counts_per_run_")

  exportTables(loaded$originalCounts,
               export.dir = paste( parameters$dir$TSV, parameters$recountID, sep = "/"),
               file.prefix = "original_counts_")

  exportTables(loaded$filtered,
               export.dir = paste( parameters$dir$TSV,parameters$recountID, sep = "/"),
               file.prefix = "filtered_counts_")


} else {
  message.with.time("Skipping load the count Table from recount experiment, merge count per sample and filter it\n","
                    from zero and near-zero variance")
}




##### Normalize the counts without log2 transformation (first test) #####
##
## Note: this method takes a table with one column per sample and one
## row per gene, we thus have to transpose the raw count table.

###### Normalization method for the recount Table after merge and filtered it ########
if (parameters$compute) {
  # dim(loaded$loaded$norm$counts)
  message.with.time("Normalizing counts based on 75th percentile")
  loaded$norm <- NormalizeCounts(
    self = loaded$filtered,
    classColumn = parameters$classColumn,
    classColors = parameters$classColor,
    # phenoTable = loaded$filteredExperiment$phenoTable,
    # classLabels = loaded$filteredExperiment$classLabels,
    method = "quantile", quantile=0.75, log2 = FALSE)
  # dim(loaded$normCounts)
  # loaded$norm$nb.samples
  # loaded$norm$nb.genes

  # class(loaded$norm)
  # loaded$norm$dataType

  # loaded$norm <- countTableWithTrainTestSets(loaded$norm)
  #  hist(unlist(loaded$loaded$norm$counts), main="Normalised count distribution", breaks=1000)

  ## Export the Normalized count tables with their associated information (pheno table, class labels) in tab-separated value (.tsv) files
  exportTables(loaded$norm,
               export.dir = paste(parameters$dir$TSV, parameters$recountID, sep = "/"),
               file.prefix = "norm_counts_")


} else {
  message.with.time("Skipping normalisation for the count Table  and log2 trasformation")
}






##### Normalize counts with log2 transformation (second test) #####
##
## Note: this method takes a table with one column per sample and one
## row per gene, we thus have to transpose the raw count table.
if (parameters$compute) {
  message.with.time("Normalizing counts based on 75th percentile + log2 transformation")
  loaded$log2norm <- NormalizeCounts(
    self = loaded$filtered,
    classColumn = parameters$classColumn,
    # counts =loaded$filteredExperiment$countTable,
    # phenoTable = loaded$filteredExperiment$phenoTable,
    # classLabels = loaded$filteredExperiment$classLabels,
    method = "quantile", quantile=0.75,
    log2 = TRUE, epsilon=0.1)

  # class(loaded$log2norm)
  print(loaded$log2norm$dataType)

  ## Export tables
  exportTables(loaded$log2norm,
               export.dir = paste(parameters$dir$TSV, parameters$recountID, sep = "/"),
               file.prefix = "log2norm_counts_")

  ## STILL IN CONSTRUCTION (2018-03-19)
  plotFigures(loaded$log2norm,
              plot.dir = file.path(dir.NormImpact),
              file.prefix = "log2norm")

}


## TO DO LATER: CHECK IF THESE FIGURES ARE WELL GENERATED, AND INCOROPORATE THEM IN THE plot.figure methods

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

#   if (ncol(loaded$log2norm$countTable) != length(loaded$log2norm$classLabels)){
#     stop(" the Number of samples in log2norm counts should be the same length of classes")
#   }
#
#
#
# } else {
#   message.with.time("Skipping normalisation for count Table with log2 transformation")
# }



## Indicate that this script has finished running
message.with.time("finished executing 02_load_and_normalise_counts.R")

##### plotting some figures to explore the nuture of recount data set #####
# message.with.time(" plotting some figures to explore distribution for the recount data set ",parameters$recountID)
# source("misc/11_impact_of_normalization_and_Log2.R")

# ##### Exhibiting the geo charactiristics for the current project #####
# message.with.time("Exhibit the geo charactiristics for such experiment: ", parameters$recountID, " in order to know the class lable
#                   for such experiment")
# head( geo.characteristics)
# geo.characteristics.file <- file.path("~/RNAseqMVA_workspace", "data", parameters$recountID, "geo.characteristics.tsv")
# write.table( geo.characteristics, file = geo.characteristics.file, quote = FALSE,
#              row.names = FALSE, sep = "\t" )
