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


  ## Define the directories where tables and figures will be stored.
  ## one directory per classifer, with separate subdirectories for tables and figures.
  classifiers <- parameters$classifiers

  #  dir.classifier <- file.path(dir.results, classifiers)
  #  names(dir.classifier) <- classifiers

  classifier.dirs <- file.path(dir.results, classifiers)
  table.dirs <- file.path(classifier.dirs, "tables")
  figure.dirs <- file.path(classifier.dirs, "figures")
  detailFigures.dir <- file.path(figure.dirs, "detailFigures")
  detailTables.dir <- file.path(table.dirs, "detailTables")
  names(detailTables.dir)<- classifiers
  names(detailFigures.dir)<- classifiers

  for (dir in c(classifier.dirs, table.dirs, figure.dirs, detailFigures.dir, detailTables.dir)) {
    #    classifier.dirs[classifier] <- file.path(dir.results, classifier)

      dir.create(dir, showWarnings = F, recursive = T)


    #    table.dirs[classifier] <- file.path(classifier.dirs[classifier], "tables")
    #dir.create(table.dirs[classifier], showWarnings = F, recursive = T)

    #    detailTables.dir[classifier] <- file.path(table.dirs[classifier], "detailTables")
    #    dir.create(detailTables.dir[classifier],showWarnings = F, recursive = T)

    #    figure.dirs[classifier] <- file.path(classifier.dirs[classifier], "figures")
    #dir.create(figure.dirs[classifier], showWarnings = F, recursive = T)

    #    detailFigures.dir[classifier] <- file.path(figure.dirs[classifier], "detailFigures")
    #    dir.create(detailFigures.dir[classifier] ,showWarnings = F, recursive = T)
  } # end loop over the dir

  ## File to store a memory image
  image.file <- file.path(dir.results, paste("RNA-seq_classifer_evaluation_", recountID, ".Rdata", sep = ""))

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


  ################################################################
  ## TO CHECK LATER: DO wE STILL NEED THESE VARIABLES ???

  ## Directory for impact of Normalization and log2 into counts (and the study of its impact)
  dir.NormImpact <- file.path(dir.results , paste("impact_of_normalisation_and_log2", sep = ""))
  dir.create(dir.NormImpact, showWarnings = F, recursive = T)

  ## Directory for the visualization of Principal component for counts (and the study of its impact)
  dir.visualisePCs <- file.path(dir.results , paste( "visualization_of_PCs", sep = ""))
  dir.create(dir.visualisePCs, showWarnings = F, recursive = T)


  # View(parameters)



  if (parameters$compute) {
    message.with.time("Loading count table from recount", "; recountID = ", parameters$recountID)
    loaded[[recountID]] <- loadCounts(recountID = recountID,
                                      mergeRuns = TRUE,
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
      initial.runs = loaded[[recountID]]$countsPerRun$nbSamples,
      initial.samples = loaded[[recountID]]$originalCounts$nbSamples,
      initial.genes = loaded[[recountID]]$originalCounts$nbGenes,
      initial.classes = loaded[[recountID]]$originalCounts$nbClasses,
      filtered.samples = loaded[[recountID]]$filtered$nbSamples,
      filtered.zerovar = length(loaded[[recountID]]$filtered$zeroVarGenes),
      filtered.nearzerovar = length(loaded[[recountID]]$filtered$nearZeroVarGenes),
      filtered.genes = loaded[[recountID]]$filtered$nbGenes,
      filtered.classes = loaded[[recountID]]$filtered$nbClasses
    )

  if (ncol(loadedStats) == 0) {
    loadedStats <- newStats
  } else {
    loadedStats <- rbind(loadedStats, newStats)
  }
}
rownames(loadedStats) <- loadedStats$recountID

## TEMPORARY: print out stats about loaded datasets
require(knitr)
kable(t(loadedStats))

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
