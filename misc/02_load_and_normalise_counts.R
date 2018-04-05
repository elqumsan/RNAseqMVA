#################### Load counts and pheno ####################
## Load a count Table from recount-experiment, merge counts per sample
## and apply some pre-filtering (remove zero-variance and near-zero-variance genes).
loaded <- list()

for (recountID in selectedRecountIDs) {

  parameters$recountID <- recountID

  ## Result directory
  dir.results <- file.path(parameters$dir$workspace, "results", parameters$recountID)


  ## Directory to exprt the tab-separate value files
  # tsv.dir <- paste(sep = "" , parameters$dir$TSV,"/",recountID)
  tsv.dir <- file.path(dir.results, "TSV")
  dir.create(path = tsv.dir, recursive = TRUE, showWarnings = FALSE)

  ## Directory to store the data downloaded from recountID
  studyPath <- file.path(parameters$dir$workspace, "data", recountID)

  ## Define the directories where tables and figures will be stored.
  ## one directory per classifer, with separate subdirectories for tables and figures.
  classifiers <- parameters$classifiers

  #  dir.classifier <- file.path(dir.results, classifiers)
  #  names(dir.classifier) <- classifiers

  classifier.dirs <- file.path(dir.results, classifiers)
  names(classifier.dirs) <- classifiers

  table.dirs <- file.path(classifier.dirs, "tables")
  names(table.dirs) <- classifiers

  figure.dirs <- file.path(classifier.dirs, "figures")
  names(figure.dirs) <- classifiers

  detailFigures.dir <- file.path(figure.dirs, "detailFigures")
  names(detailFigures.dir) <- classifiers

  detailTables.dir <- file.path(table.dirs, "detailTables")
  names(detailTables.dir) <- classifiers

  ## Create all the recountID-specific sub-directories
  for (dir in c(classifier.dirs, table.dirs, figure.dirs, detailFigures.dir, detailTables.dir)) {
      dir.create(dir, showWarnings = F, recursive = T)
  } # end loop over the dir

  ################################################################
  ## TO CHECK LATER: DO wE STILL NEED THESE VARIABLES ???

  ## Directory for impact of Normalization and log2 into counts (and the study of its impact)
  dir.NormImpact <- file.path(dir.results , paste("impact_of_normalisation_and_log2", sep = ""))
  dir.create(dir.NormImpact, showWarnings = F, recursive = T)

  ## Directory for the visualization of Principal component for counts (and the study of its impact)
  dir.visualisePCs <- file.path(dir.results , paste( "visualization_of_PCs", sep = ""))
  dir.create(dir.visualisePCs, showWarnings = F, recursive = T)


  if (parameters$reload == TRUE) {
    ################################################################################
    ## Save an image of the memory, so I can reload it later to avoid re-running all the analyses
    parameters.current <- parameters # Keep current parameters to restore them after having loaded a memory image
    message.with.time("Loading memory image ")
    load(file = image.file)
    parameters <- parameters.current ## Reload current parameters (they might have been saved different in the memory image)
    rm(parameters.current)
  }



  ## Overwrite default parameters wih project-specific parameters
  selected.parameters <- project.parameters[[recountID]]
  if (is.null(selected.parameters)) {
    message("No specific parameters for recount ID ", recountID)
    message("Using generic parameters from the yaml file. ")
  } else {
    message("Using specific parameters specfied in yaml file for recount ID ", recountID)
    parameters[names(selected.parameters)] <- project.parameters[[recountID]]
    names(parameters$data.types)<-parameters$data.types
    names(parameters$variables.type)<-parameters$variables.type
  }

  ## Prefix for experiments with permuted class labels
  ## (negative controls to estimate random expectation)
  perm.prefix <- parameters$perm_prefix



  # View(parameters)



  if (parameters$compute) {
    message.with.time("Loading count table from recount", "; recountID = ", parameters$recountID)
    loaded[[recountID]] <- loadCounts(recountID = recountID,
                                      mergeRuns = parameters$mergeRuns,
                                      classColumn = parameters$classColumn,
                                      minSamplesPerClass = parameters$minSamplesPerClass,
                                      na.rm = parameters$na.rm)


    ## Select training and testing sets on the filtered table with raw counts
    ## These wil then be passed to all the derived count tables (normalised, DGE, ...)
    loaded[[recountID]]$filtered <- countTableWithTrainTestSets(loaded[[recountID]]$filtered)

    #### Export the count tables with their associated information (pheno table, class labels) in tab-separated value (.tsv) files ###
    exportTables(loaded[[recountID]]$countsPerRun,
                 export.dir = paste(parameters$dir$TSV, parameters$recountID, sep = "/"),
                 file.prefix = "counts_per_run_")

    exportTables(loaded[[recountID]]$originalCounts,
                 export.dir = paste( parameters$dir$TSV, parameters$recountID, sep = "/"),
                 file.prefix = "original_counts_")

    exportTables(loaded[[recountID]]$filtered,
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
    # dim(loaded[[recountID]]$loaded$norm$counts)
    message.with.time("Normalizing counts based on 75th percentile")
    loaded[[recountID]]$norm <- NormalizeCounts(
      self = loaded[[recountID]]$filtered,
      classColumn = parameters$classColumn,
      classColors = parameters$classColor,
      # phenoTable = loaded[[recountID]]$filteredExperiment$phenoTable,
      # classLabels = loaded[[recountID]]$filteredExperiment$classLabels,
      method = "quantile", quantile=0.75, log2 = FALSE)
    # dim(loaded[[recountID]]$normCounts)
    # loaded[[recountID]]$norm$nb.samples
    # loaded[[recountID]]$norm$nb.genes

    # class(loaded[[recountID]]$norm)
    # loaded[[recountID]]$norm$dataType

    # loaded[[recountID]]$norm <- countTableWithTrainTestSets(loaded[[recountID]]$norm)
    #  hist(unlist(loaded[[recountID]]$loaded[[recountID]]$norm$counts), main="Normalised count distribution", breaks=1000)

    ## Export the Normalized count tables with their associated information (pheno table, class labels) in tab-separated value (.tsv) files
    exportTables(loaded[[recountID]]$norm,
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
    loaded[[recountID]]$log2norm <- NormalizeCounts(
      self = loaded[[recountID]]$filtered,
      classColumn = parameters$classColumn,
      # counts =loaded$filteredExperiment$countTable,
      # phenoTable = loaded$filteredExperiment$phenoTable,
      # classLabels = loaded$filteredExperiment$classLabels,
      method = "quantile", quantile=0.75,
      log2 = TRUE, epsilon=0.1)

    # class(loaded[[recountID]]$log2norm)
    # print(loaded[[recountID]]$log2norm$dataType)

    ## Export tables
    exportTables(loaded[[recountID]]$log2norm,
                 export.dir = paste(parameters$dir$TSV, parameters$recountID, sep = "/"),
                 file.prefix = "log2norm_counts_")

    # ## STILL IN CONSTRUCTION (2018-03-19)

    # plotFigures(loaded[[recountID]]$log2norm,

    # plotFigures(loaded$log2norm,

    #             plot.dir = file.path(dir.NormImpact),
    #             file.prefix = "log2norm")

  }


  ## TO DO LATER: CHECK IF THESE FIGURES ARE WELL GENERATED, AND INCOROPORATE THEM IN THE plot.figure methods

  # plot.file <- file.path(dir.NormImpact, "log2normCount_hist.pdf")
  # message("\tlog2(norm counts) histogram\t", plot.file)
  # pdf(plot.file, width=7, height=5)
  # hist(unlist(loaded[[recountID]]$log2norm$counts), breaks=100,
  #      col="grey",
  #      main=paste("log2(norm counts) distrib;", recountID),
  #      las=1,
  #      xlab="log2(norm counts)",
  #      ylab="Frequency")
  #
  # silence <- dev.off()

  #   if (ncol(loaded[[recountID]]$log2norm$countTable) != length(loaded[[recountID]]$log2norm$classLabels)){
  #     stop(" the Number of samples in log2norm counts should be the same length of classes")
  #   }
  #
  #
  #
  # } else {
  #   message.with.time("Skipping normalisation for count Table with log2 transformation")
  # }


  # Check loaded[[recountID]] objects

  # attributes(loaded[[recountID]])
  # class(loaded[[recountID]]$countsPerRun)
  # class(loaded[[recountID]]$originalCounts)
  # class(loaded[[recountID]]$filtered)
  # class(loaded[[recountID]]$norm)
  # class(loaded[[recountID]]$log2norm)
  #
  #
  # # unlist(lapply(loaded[[recountID]]$norm$trainTestProperties$trainIndices, length))
  # length(unlist(loaded[[recountID]]$norm$trainTestProperties$trainIndices))
  # sum(unlist(loaded[[recountID]]$norm$trainTestProperties$trainIndices) != unlist(loaded[[recountID]]$filtered$trainTestProperties$trainIndices))
  # sum(unlist(loaded[[recountID]]$norm$trainTestProperties$trainIndices) != unlist(loaded[[recountID]]$log2norm$trainTestProperties$trainIndices))


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

} # end loop over recountIDs


## Compute statistics about loaded datasets
loadedStats <- data.frame()
for (recountID in selectedRecountIDs) {
  newStats <-
    data.frame(
      recountID = recountID,
      runs = loaded[[recountID]]$countsPerRun$nbSamples,
      samples = loaded[[recountID]]$originalCounts$nbSamples,
      samples.filtered = loaded[[recountID]]$filtered$nbSamples,

      genes.ori = loaded[[recountID]]$originalCounts$nbGenes,
      genes.NA = length(loaded[[recountID]]$filtered$naGenes),
      genes.zeroVar = length(loaded[[recountID]]$filtered$zeroVarGenes),
      genes.NZfilter = length(loaded[[recountID]]$filtered$nearZeroVarGenes),
      genes.filtered = loaded[[recountID]]$filtered$nbGenes,

      classes.ori = loaded[[recountID]]$originalCounts$nbClasses,
      classes.filtered = loaded[[recountID]]$filtered$nbClasses
    )

  if (ncol(loadedStats) == 0) {
    loadedStats <- newStats
  } else {
    loadedStats <- rbind(loadedStats, newStats)
  }
}
rownames(loadedStats) <- loadedStats$recountID


#loadedStats$genes.NA + loadedStats$genes.zeroVar + loadedStats$genes.NZfilter + loadedStats$genes.filtered
#loadedStats$genes.ori

loadedStats$pc.NA <-100*loadedStats$genes.NA / loadedStats$genes.ori
loadedStats$pc.zeroVar <- 100*loadedStats$genes.zeroVar / loadedStats$genes.ori
loadedStats$pc.NZfilter <- 100*loadedStats$genes.NZfilter / loadedStats$genes.ori
loadedStats$pc.kept <-  100*loadedStats$genes.filtered / loadedStats$genes.ori

## TEMPORARY: print out stats about loaded datasets
require(knitr)
kable(t(loadedStats))

write.table(x = t(loadedStats), file = file.path(parameters$dir$results, "experiment_summaries.tsv"),
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)

## YOUR MISSION: GENERATE A BARPLOT or a set of piecharts SHOWING THE FOLLOWING NUMBERS for the different study cases
## - NA(to be added to the result of filtering)
## - zero var
## - near zero var
## - kept genes
## The total shoud give the same result as initial.genes

gene.pc <- loadedStats[, c("pc.NA", "pc.zeroVar", "pc.NZfilter", "pc.kept")]
row.names(gene.pc) <- loadedStats$recountID
## Order experiments by increasing percentof kept genes (so that the highest come on top of the barplot)
gene.pc <- gene.pc[order(gene.pc$pc.kept, decreasing = FALSE), ]
## apply(gene.pc, 1, sum) ## The sums must give 100 for each experiment

file.prefix <- paste("experiments_sumaries.pdf")
barPlot.file <- file.path(parameters$dir$results,file.prefix)
message("Filtering summary barplot: ", barPlot.file)
pdf(file = barPlot.file, width=7, height=2+1*nrow(loadedStats))
save.margins <- par("mar")
par(mar= c(5,7,5,1))

#kept.label <- paste(sep="",round(digits=0, gene.pc$pc.kept), "%")
ypos <- barplot(t(gene.pc),las=1, horiz = TRUE,
                 col = c("black", "red", "orange", "#44DD44"),
                 legend.text = c("NA values", "Zero var", "NZ filter", "kept"),
                 main = "Filtering impact on study cases",
                 xlab = "Proportions of genes",
#                 names.arg = paste(sep="", rownames(gene.pc), " (", kept.label, ")"),
                 xlim=c(0, 170))
text(x = 100, kept.label, y = ypos, pos = 2)

par(mar = save.margins)
silence<- dev.off()

#### Draw a barplot with the number of samples per class ####
file.prefix <- paste("samples_per_classes.pdf")
barPlot.file <- file.path(parameters$dir$results,file.prefix)
message("Filtering summary barplot: ", barPlot.file)
pdf(file = barPlot.file, width=8, height=12)
save.margins <- par("mar")
par(mfrow=c(4,2))
par(mar=c(4,15,5,1))
for (recountID in names(loaded)) {
  heights <- barplot(sort(loaded[[recountID]]$original$samplesPerClass, decreasing = TRUE),
          horiz = TRUE, las=1, cex.names = 0.7, main=recountID, xlab="Samples per class")
  barplot(sort(loaded[[recountID]]$filtered$samplesPerClass, decreasing = TRUE), add=TRUE,
          horiz = TRUE, las=1, cex.names = 0.7, main=recountID, xlab="Samples per class", col="#00BB00")
}
par(mfrow=c(1,1))
par(mar = save.margins)
silence<- dev.off()

## Save a memory image that can be re-loaded next time to avoid re-computing all the normalisation and so on.
if (parameters$save.image) {
  dir.create(parameters$dir$memoryImages, showWarnings = FALSE, recursive = TRUE)
  loaded.mem.image <- file.path(parameters$dir$memoryImages, "data_loaded.Rdata")
  message.with.time("Saving memory image after loading: ", loaded.mem.image)
  save.image(file = loaded.mem.image)
}

## Indicate that this script has finished running
message.with.time("finished executing 02_load_and_normalise_counts.R")

#######################################################
# loadedObjects.dataType <- data.frame()
#
# for(recountID in selectedRecountIDs){
#   loadedDataType <-
#     data.frame(
#       "originalCounts" = loaded[[recountID]]$originalCounts ,
#       "filtere" = loaded[[recountID]]$filtered,
#       "norm" = loaded[[recountID]]$norm,
#       "log2norm" = loaded[[recountID]]$log2norm )
#
#   if (ncol(loadedDataType)== 0){
#     loadedObjects.dataType <- loadedDataType
#   }else {
#     loadedObjects.dataType <-rbind(loadedObjects.dataType, loadedDataType)
#
#   }
#
# }
#
