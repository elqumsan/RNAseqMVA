#################### Load counts and pheno ####################
## Load a count Table from recount-experiment, merge counts per sample
## and apply some pre-filtering (remove zero-variance and near-zero-variance genes).

studyCases <- list() ## a list containing all the loaded datasets + their pre-processed data

for (recountID in selectedRecountIDs) {

  #### Specify generic and recountID-specific parameters ####
  parameters <- initRecountID(recountID, configFile)

  # ## Load default parameters for each new recountID
  # ## (was previously parsed from the YAML file)
  # parameters <- project.parameters$default
  #
  # ## Specify the current recountID in parameters
  # parameters$recountID <- recountID
  #
  # ## Overwrite default parameters wih project-specific parameters
  # selected.parameters <- project.parameters[[recountID]]
  # if (is.null(selected.parameters)) {
  #   message("No specific parameters for recount ID ", recountID)
  #   message("Using generic parameters from the yaml file. ")
  # } else {
  #   message("Using specific parameters specfied in yaml file for recount ID ", recountID)
  #   parameters[names(selected.parameters)] <- project.parameters[[recountID]]
  #   names(parameters$data.types)<-parameters$data.types
  #   names(parameters$variables.type)<-parameters$variables.type
  # }
  #
  # ## Convert list-formatted class colors to named vector (YAML does not allow to specify named vectors)
  # if (!is.null(parameters$classColors)) {
  #   if (class(parameters$classColors) == "list") {
  #     parameters$classColors <- unlist(parameters$classColors)
  #   }
  # }

  ## BEWARE: NOT SURE THIS IS FUNCTIONAL, RELOADING MEMORY IMAGES WILL BE TREATED LATER
  # if (parameters$reload == TRUE) {
  #   ## Save an image of the memory, so I can reload it later to avoid re-running all the analyses
  #   parameters.current <- parameters # Keep current parameters to restore them after having loaded a memory image
  #   message.with.time("Loading memory image ")
  #   load(file = image.file)
  #   parameters <- parameters.current ## Reload current parameters (they might have been saved different in the memory image)
  #   rm(parameters.current)
  # }


  ## Prefix for experiments with permuted class labels
  ## (negative controls to estimate random expectation)
  perm.prefix <- parameters$perm_prefix


  #### Specification of the directories   ####

  ## MUSTAFA: TO DO SOME DAY: THE DIRECTORIES SHOULD BE ATTRIBUTES OF THE OBJECTS RATHER THAN GLOBAL VARIABLES

  # Main directory should be adapted to the user's configuration
  dir.main <- parameters$dir$main

  ## All other directories should be defined relative to dir.main
  # dir.scripts <- file.path(dir.main, "R")

  # ## Result directory
  # parameters$dir$results <- file.path(parameters$dir$workspace, "results", parameters$recountID)
  #
  ## Directory to exprt the tab-separate value files
  # parameters$dir$tsv <- paste(sep = "" , parameters$dir$tsv,"/",recountID)
  # parameters$dir$tsv <- file.path(parameters$dir$results, "TSV")
  # dir.create(path = parameters$dir$tsv, recursive = TRUE, showWarnings = FALSE)

  ## Directory to store the data downloaded from recountID
  studyPath <- file.path(parameters$dir$workspace, "data", recountID)



  # View(parameters)



  if (parameters$compute) {
    message.with.time("Loading count table from recount", "; recountID = ", parameters$recountID)
    studyCases[[recountID]] <- loadCounts(recountID = recountID,
                                      mergeRuns = parameters$mergeRuns,
                                      classColumn = parameters$classColumn,
                                      minSamplesPerClass = parameters$minSamplesPerClass,
                                      na.rm = parameters$na.rm)

    ## Attach the recountID-specific parameters to the loaded data.
    studyCases[[recountID]]$parameters <- parameters

    ## Select training and testing sets on the filtered table with raw counts
    ## These wil then be passed to all the derived count tables (normalised, DGE, ...)
    studyCases[[recountID]]$filtered <- countTableWithTrainTestSets(studyCases[[recountID]]$filtered)

    #### Export the count tables with their associated information (pheno table, class labels) in tab-separated value (.tsv) files ###
    exportTables(studyCases[[recountID]]$countsPerRun,
                 export.dir = file.path(parameters$dir$tsv, parameters$recountID),
                 file.prefix = "counts_per_run_")

    exportTables(studyCases[[recountID]]$originalCounts,
                 export.dir = file.path(parameters$dir$tsv, parameters$recountID),
                 file.prefix = "original_counts_")

    exportTables(studyCases[[recountID]]$filtered,
                 export.dir = file.path(parameters$dir$tsv,parameters$recountID),
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
    # dim(studyCases[[recountID]]$studyCases$norm$counts)
    message.with.time("Normalizing counts based on 75th percentile")
    studyCases[[recountID]]$norm <- NormalizeCounts(
      self = studyCases[[recountID]]$filtered,
      classColumn = parameters$classColumn,
      classColors = parameters$classColor,
      # phenoTable = studyCases[[recountID]]$filteredExperiment$phenoTable,
      # classLabels = studyCases[[recountID]]$filteredExperiment$classLabels,
      method = "quantile", quantile=0.75, log2 = FALSE)
    # dim(studyCases[[recountID]]$normCounts)
    # studyCases[[recountID]]$norm$nb.samples
    # studyCases[[recountID]]$norm$nb.genes

    # class(studyCases[[recountID]]$norm)
    # studyCases[[recountID]]$norm$dataType

    # studyCases[[recountID]]$norm <- countTableWithTrainTestSets(studyCases[[recountID]]$norm)
    #  hist(unlist(studyCases[[recountID]]$studyCases[[recountID]]$norm$counts), main="Normalised count distribution", breaks=1000)

    ## Export the Normalized count tables with their associated information (pheno table, class labels) in tab-separated value (.tsv) files
    exportTables(studyCases[[recountID]]$norm,
                 export.dir = paste(parameters$dir$tsv, parameters$recountID, sep = "/"),
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
    studyCases[[recountID]]$log2norm <- NormalizeCounts(
      self = studyCases[[recountID]]$filtered,
      classColumn = parameters$classColumn,
      # counts =studyCases$filteredExperiment$countTable,
      # phenoTable = studyCases$filteredExperiment$phenoTable,
      # classLabels = studyCases$filteredExperiment$classLabels,
      method = "quantile", quantile=0.75,
      log2 = TRUE, epsilon=0.1)

    # class(studyCases[[recountID]]$log2norm)
    # print(studyCases[[recountID]]$log2norm$dataType)

    ## Export tables
    exportTables(studyCases[[recountID]]$log2norm,
                 export.dir = paste(parameters$dir$tsv, parameters$recountID, sep = "/"),
                 file.prefix = "log2norm_counts_")

    # ## STILL IN CONSTRUCTION (2018-03-19)

    # plotFigures(studyCases[[recountID]]$log2norm,

    # plotFigures(studyCases$log2norm,

    #             plot.dir = file.path(parameters$dir$NormalizationImpact),
    #             file.prefix = "log2norm")

  }


  ## TO DO LATER: CHECK IF THESE FIGURES ARE WELL GENERATED, AND INCOROPORATE THEM IN THE plot.figure methods

  # plot.file <- file.path(parameters$dir$NormalizationImpact, "log2normCount_hist.pdf")
  # message("\tlog2(norm counts) histogram\t", plot.file)
  # pdf(plot.file, width=7, height=5)
  # hist(unlist(studyCases[[recountID]]$log2norm$counts), breaks=100,
  #      col="grey",
  #      main=paste("log2(norm counts) distrib;", recountID),
  #      las=1,
  #      xlab="log2(norm counts)",
  #      ylab="Frequency")
  #
  # silence <- dev.off()

  #   if (ncol(studyCases[[recountID]]$log2norm$countTable) != length(studyCases[[recountID]]$log2norm$classLabels)){
  #     stop(" the Number of samples in log2norm counts should be the same length of classes")
  #   }
  #
  #
  #
  # } else {
  #   message.with.time("Skipping normalisation for count Table with log2 transformation")
  # }


  # Check studyCases[[recountID]] objects

  # attributes(studyCases[[recountID]])
  # class(studyCases[[recountID]]$countsPerRun)
  # class(studyCases[[recountID]]$originalCounts)
  # class(studyCases[[recountID]]$filtered)
  # class(studyCases[[recountID]]$norm)
  # class(studyCases[[recountID]]$log2norm)
  #
  #
  # # unlist(lapply(studyCases[[recountID]]$norm$trainTestProperties$trainIndices, length))
  # length(unlist(studyCases[[recountID]]$norm$trainTestProperties$trainIndices))
  # sum(unlist(studyCases[[recountID]]$norm$trainTestProperties$trainIndices) != unlist(studyCases[[recountID]]$filtered$trainTestProperties$trainIndices))
  # sum(unlist(studyCases[[recountID]]$norm$trainTestProperties$trainIndices) != unlist(studyCases[[recountID]]$log2norm$trainTestProperties$trainIndices))


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
studyCasesStats <- data.frame()

for (recountID in selectedRecountIDs) {
  newStats <-
    data.frame(
      recountID = recountID,
      runs = studyCases[[recountID]]$countsPerRun$nbSamples,
      samples = studyCases[[recountID]]$originalCounts$nbSamples,
      samples.filtered = studyCases[[recountID]]$filtered$nbSamples,

      genes.ori = studyCases[[recountID]]$originalCounts$nbGenes,
      genes.NA = length(studyCases[[recountID]]$filtered$naGenes),
      genes.zeroVar = length(studyCases[[recountID]]$filtered$zeroVarGenes),
      genes.NZfilter = length(studyCases[[recountID]]$filtered$nearZeroVarGenes),
      genes.filtered = studyCases[[recountID]]$filtered$nbGenes,

      classes.ori = studyCases[[recountID]]$originalCounts$nbClasses,
      classes.filtered = studyCases[[recountID]]$filtered$nbClasses
    )

  if (ncol(studyCasesStats) == 0) {
    studyCasesStats <- newStats
  } else {
    studyCasesStats <- rbind(studyCasesStats, newStats)
  }
}
rownames(studyCasesStats) <- studyCasesStats$recountID


#studyCasesStats$genes.NA + studyCasesStats$genes.zeroVar + studyCasesStats$genes.NZfilter + studyCasesStats$genes.filtered
#studyCasesStats$genes.ori

studyCasesStats$pc.NA <-100*studyCasesStats$genes.NA / studyCasesStats$genes.ori
studyCasesStats$pc.zeroVar <- 100*studyCasesStats$genes.zeroVar / studyCasesStats$genes.ori
studyCasesStats$pc.NZfilter <- 100*studyCasesStats$genes.NZfilter / studyCasesStats$genes.ori
studyCasesStats$pc.kept <-  100*studyCasesStats$genes.filtered / studyCasesStats$genes.ori

## TEMPORARY: print out stats about loaded datasets
require(knitr)
kable(t(studyCasesStats))

write.table(x = t(studyCasesStats), file = file.path(parameters$dir$results, "experiment_summaries.tsv"),
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)

## YOUR MISSION: GENERATE A BARPLOT or a set of piecharts SHOWING THE FOLLOWING NUMBERS for the different study cases
## - NA(to be added to the result of filtering)
## - zero var
## - near zero var
## - kept genes
## The total shoud give the same result as initial.genes

gene.pc <- studyCasesStats[, c("pc.NA", "pc.zeroVar", "pc.NZfilter", "pc.kept")]
row.names(gene.pc) <- studyCasesStats$recountID
## Order experiments by increasing percentof kept genes (so that the highest come on top of the barplot)
gene.pc <- gene.pc[order(gene.pc$pc.kept, decreasing = FALSE), ]
## apply(gene.pc, 1, sum) ## The sums must give 100 for each experiment

figure.file <- paste("experiments_sumaries.pdf")
barPlot.file <- file.path(parameters$dir$results,figure.file)
message("Filtering summary barplot: ", barPlot.file)
pdf(file = barPlot.file, width=7, height=2+1*nrow(studyCasesStats))
save.margins <- par("mar")
par(mar= c(5,7,5,1))

kept.label <- paste(sep="",round(digits=0, gene.pc$pc.kept), "%")
ypos <- barplot(t(gene.pc), las=1, horiz = TRUE,
                col = c("black", "red", "orange", "#44DD44"),
                legend.text = c("NA values", "Zero var", "NZ filter", "kept"),
                main = "Filtering impact on study cases",
                xlab = "Proportions of genes",
                #                 names.arg = paste(sep="", rownames(gene.pc), " (", kept.label, ")"),
                xlim=c(0, 170))
text(x = 100, kept.label, y = ypos, pos = 4)

par(mar = save.margins)
silence<- dev.off()

#### Draw a barplot with the number of samples per class ####
figure.file <- paste("samples_per_classes.pdf")
barPlot.file <- file.path(parameters$dir$results,figure.file)
message("Filtering summary barplot: ", barPlot.file)
pdf(file = barPlot.file, width=8, height=12)
save.margins <- par("mar")
par(mfrow=c(4,2))
par(mar=c(4,15,5,1))
for (recountID in names(studyCases)) {
  heights <- barplot(
    sort(studyCases[[recountID]]$original$samplesPerClass, decreasing = TRUE),
    horiz = TRUE, las=1, cex.names = 0.7, main=recountID,
    xlab="Samples per class", col="white")
  barplot(sort(studyCases[[recountID]]$filtered$samplesPerClass, decreasing = TRUE),
          add=TRUE, horiz = TRUE, las=1, cex.names = 0.7,
          main=recountID, xlab="Samples per class",
          col = studyCases[[recountID]]$filtered$classColors)
#          col="#00BB00")
  abline(v = studyCases[[recountID]]$parameters$minSamplesPerClass, col="red")
}
par(mfrow=c(1,1))
par(mar = save.margins)
silence<- dev.off()

## Save a memory image that can be re-loaded next time to avoid re-computing all the normalisation and so on.
if (parameters$save.image) {
  dir.create(parameters$dir$memoryImages, showWarnings = FALSE, recursive = TRUE)
  studyCases.mem.image <- file.path(parameters$dir$memoryImages, "data_studyCases.Rdata")
  message.with.time("Saving memory image after loading: ", studyCases.mem.image)
  save.image(file = studyCases.mem.image)
}

## Indicate that this script has finished running
message.with.time("finished executing 02_load_and_normalise_counts.R")

#######################################################
# loadedObjects.dataType <- data.frame()
#
# for(recountID in selectedRecountIDs){
#   loadedDataType <-
#     data.frame(
#       "originalCounts" = studyCases[[recountID]]$originalCounts ,
#       "filtere" = studyCases[[recountID]]$filtered,
#       "norm" = studyCases[[recountID]]$norm,
#       "log2norm" = studyCases[[recountID]]$log2norm )
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

