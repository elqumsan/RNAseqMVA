#################### Load counts and pheno ####################
## Load a count Table from recount-experiment, merge counts per sample
## and apply some pre-filtering (remove zero-variance and near-zero-variance genes).

message("Loading study cases")


#### Check the recountID ####
studyCases <- list() ## a list containing all the loaded datasets + their pre-processed data

recountID <- selectedRecountIDs[1] ## We only load one recountID at a time

if (length(selectedRecountIDs) > 1) {
  stop("Multiple study case analysis is not supported, please specify a single study case")
}

#### Init the generic and recountID-specific parameters ####
parameters <- initRecountID(recountID, project.parameters)
# View(parameters)

message.with.time("Building StudyCase for recountID\t", recountID,
                  "\t", parameters$short_label,
                  "\n\tfeature type: ", project.parameters$global$feature)

#### Load study case ####
studyCases[[recountID]] <- StudyCase(recountID = recountID, parameters = parameters)
message("Finished loading study case ", recountID,
        "\t", parameters$short_label,
        "\n\tfeature type: ", parameters$feature)


#### Export the count tables with their associated information (pheno table, class labels) in tab-separated value (.tsv) files #####
if (project.parameters$global$export.tables) {
  exportTables(studyCases[[recountID]])
}
message("Finished exporting data table")


#### Plot histograms of log2 normalized counts ####
datasetNames <- names(studyCases[[recountID]]$datasetsForTest)
log2countNames <- grep(pattern = "log2$", datasetNames, value = TRUE)
shortLabel <- parameters$short_label
for (datasetName in log2countNames) {
  plot.file <- file.path(parameters$dir$NormalizationImpact,
                         paste0(recountID, "_", datasetName, "_hist", ".pdf"))
  message("\t", datasetName, " histogram", "\t", plot.file)
  pdf(plot.file, width = 7, height = 5)
  par.ori <- par(no.readonly = TRUE)
  par(mar = c(5.1, 6.1, 4.1, 1))
  hist(unlist(studyCases[[recountID]]$datasetsForTest[[datasetName]]$dataTable), breaks = 100,
       col = "grey",
       main = paste0(datasetName, " count distribution",
                     "\n", shortLabel, "; ", recountID),
       las = 1,
       xlab = "log2(norm counts)",
       ylab = NA)
  par <- par.ori
  silence <- dev.off()
}

#### Save the study case object in Rdata format ####
## This is useful to save the study case in a memory image that can be re-loaded later
## Reload can be done with the script misc/02b_reload_counts.R
if (project.parameters$global$save.image) {
  studyCase <- studyCases[[recountID]]
  featureType <- studyCase$parameters$feature
  studyCaseFile <-  file.path(
    project.parameters$global$dir$memoryImages,
    paste0(
      recountID,
      "_", featureType,
      "_loaded_studyCase",
      ".Rdata"))

  message("Saving study case\t", recountID, "\t", studyCaseFile)
  save(studyCase, file = studyCaseFile)

  for (datasetName in names(studyCase$datasetsForTest)) {
    dataset <- studyCase$datasetsForTest[[datasetName]]
    featureType <- dataset$parameters$feature
    datasetFile <-  file.path(
      project.parameters$global$dir$memoryImages,
      paste0(recountID, "_", featureType, "_loaded-dataset_", datasetName, ".Rdata"))
    message("Saving dataset\t", recountID, "\t", datasetName, "\t", datasetFile)
    save(dataset, file = datasetFile)

  }
  rm(dataset)
  rm(studyCase)
}

message("Finished exporting memory image of studyCase object ", length(studyCases))



#### Compute statistics about loaded datasets ####
studyCasesStats <- data.frame()
message("Computing summary statistics for ", length(studyCases), " study case(s). ")
for (recountID in selectedRecountIDs) {
  newStats <-
    data.frame(
      recountID = recountID,
      runs = studyCases[[recountID]]$rawData$countsPerRun$nbSamples,
      samples = studyCases[[recountID]]$rawData$countsPerSample$nbSamples,
      samples.filtered = studyCases[[recountID]]$datasetsForTest$filtered$nbSamples,
      classes.ori = studyCases[[recountID]]$rawData$countsPerSample$nbClasses,
      classes.filtered = studyCases[[recountID]]$datasetsForTest$filtered$nbClasses
    )


  featureType <- parameters$feature
  newStats[paste0(featureType, ".ori")] <- studyCases[[recountID]]$rawData$countsPerSample$nbGenes
  newStats[paste0(featureType, ".NA")] <- length(studyCases[[recountID]]$datasetsForTest$filtered$naGenes)
  newStats[paste0(featureType, ".zeroVar")] <- length(studyCases[[recountID]]$datasetsForTest$filtered$zeroVarGenes)
  newStats[paste0(featureType, ".NZfilter")] <- length(studyCases[[recountID]]$datasetsForTest$filtered$nearZeroVarGenes)
  newStats[paste0(featureType, ".filtered")] <- studyCases[[recountID]]$datasetsForTest$filtered$nbGenes


  if (ncol(studyCasesStats) == 0) {
    studyCasesStats <- newStats
  } else {
    studyCasesStats <- rbind(studyCasesStats, newStats)
  }
}
rownames(studyCasesStats) <- studyCasesStats$recountID
studyCasesStats$pc.NA <- 100*studyCasesStats[paste0(featureType, ".NA")] / studyCasesStats[paste0(featureType, ".ori")]
studyCasesStats$pc.zeroVar <- 100*studyCasesStats[paste0(featureType, ".zeroVar")] / studyCasesStats[paste0(featureType, ".ori")]
studyCasesStats$pc.NZfilter <- 100*studyCasesStats[paste0(featureType, ".NZfilter")] / studyCasesStats[paste0(featureType, ".ori")]
studyCasesStats$pc.kept <-  100*studyCasesStats[paste0(featureType, ".filtered")] / studyCasesStats[paste0(featureType, ".ori")]

## DEBUG: print out stats about loaded datasets
# library(knitr)
# kable(t(as.data.frame(studyCasesStats)))

## Export the studyCase statistics in a TSV file
studyCaseStatsFile <- file.path(
  parameters$dir$tsv,
  paste0(recountID, "_", parameters$feature, "_experiment_summaries.tsv"))
write.table(x = t(studyCasesStats),
            file = studyCaseStatsFile,
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)
message("\tstudy case statistics saved in file\t", studyCaseStatsFile)


####  Feature filtering barplots:  proportions of features filtered out or kept in each study case ####
## NOTE: this plot was generated at the time when we were loading all the study cases together,
## which we don't do anymore.
## This is not tractable (memory) so we comment it

# gene.pc <- studyCasesStats[, c("pc.NA", "pc.zeroVar", "pc.NZfilter", "pc.kept")]
# row.names(gene.pc) <- studyCasesStats$recountID

## Order experiments by increasing percentof kept genes (so that the highest come on top of the barplot)
# gene.pc <- gene.pc[order(gene.pc$pc.kept, decreasing = FALSE), ]
## apply(gene.pc, 1, sum) ## The sums must give 100 for each experiment
# figure.file <- paste0(featureType, "_filtering_summaries.pdf")
# barPlot.file <- file.path(parameters$dir$figures, figure.file)
# message("Filtering summary barplot: ", barPlot.file)
# pdf(file = barPlot.file, width = 7, height = 2 + 1 * nrow(studyCasesStats))
# save.margins <- par("mar")
# par(mar = c(5,7,5,1))

# kept.label <- paste0(round(digits = 0, gene.pc$pc.kept), "%")
# ypos <- barplot(t(gene.pc), las = 1, horiz = TRUE,
#                 col = c("black", "red", "orange", "#44DD44"),
#                 legend.text = c("NA values", "Zero var", "NZ filter", "kept"),
#                 main = "Filtering impact on study cases",
#                 xlab = "Proportions of genes",
#                 xlim = c(0, 170))
# text(x = 100, kept.label, y = ypos, pos = 4)

# par(mar = save.margins)  ## Restore original graphical parameter
# silence <- dev.off()      ## Close graphical device

#### Draw a barplot with the number of samples per class ####
message("Exporting barplots with number of samples per class")
for (recountID in names(studyCases)) {
#  dir.figures <- file.path(parameters$dir$results, recuntID, "figures")
#  dir.create(dir.figures, recursive = TRUE, showWarnings = FALSE)
  figure.file <- paste0(recountID, "_", featureType, "_samples_per_classes.pdf")
  barPlot.file <- file.path(parameters$dir$figures, figure.file)
  message("\t", recountID, "\tSamples per classes barplot: ", barPlot.file)
  pdf(file = barPlot.file, width = 5, height = 8)
  par.ori <- par(no.readonly = TRUE)
  # par(mfrow = c(4,2))
  classNames <- studyCases[[recountID]]$rawData$countsPerSample$classNames
  classNameLen <- max(nchar(classNames))
  par(mar = c(4, 1 + classNameLen*0.3, 5, 1))
  # par("mar")
  # heights <- barplot(
  #   sort(studyCases[[recountID]]$rawData$countsPerSample$samplesPerClass, decreasing = TRUE),
  #   horiz = TRUE, las = 1, cex.names = 0.7, main = recountID,
  #   xlab = "Samples per class", col = "white")
  barplot(sort(studyCases[[recountID]]$datasetsForTest$filtered$samplesPerClass, decreasing = TRUE),
#          add = TRUE,
          horiz = TRUE, las = 1, cex.names = 0.7,
#          main = "BOUM",
          main = paste0(recountID, " ", featureType, "\n", parameters$short_label),
          xlab = "Samples per class",
          col = studyCases[[recountID]]$datasetsForTest$filtered$classColors)
  abline(v = studyCases[[recountID]]$parameters$filtering$minSamplesPerClass, col = "red")
  par(mfrow = c(1,1))
  par(par.ori)
  silence <- dev.off()
}


#### Principal component plots ####
message("\n\tExporting principal component plots")
for (recountID in selectedRecountIDs) {
  studyCase <- studyCases[[recountID]]
  parameters <- studyCase$parameters
  featureType <- parameters$feature
  message("\tPC plots for study case ", recountID,
          "\t", featureType,
          "\t", studyCase$parameters$short_label )

  plotDir <- parameters$dir$PCviz
  dir.create(plotDir, recursive = TRUE, showWarnings = FALSE)
  message("\t\tPC plots saved in directory\t", plotDir)

  ## Identify datasets corresponding to PCs
  datasetNames <- names(studyCase$datasetsForTest)
  pcNames <- grep(datasetNames, pattern = "_PC$", perl = TRUE, value = TRUE)

  ## Variance barplot of the components
  # datasetName <- "q0.75_log2_PC"
  for (datasetName in pcNames) {
    dataset <- studyCase$datasetsForTest[[datasetName]]

    PCplot.file <- file.path(plotDir, paste0(recountID, "_", featureType, "_", datasetName, "_variance.pdf"))
    message("\t\tPC variance plot: ", PCplot.file)
    pdf(file = PCplot.file, width = 7, height = 5)
    plot(dataset$prcomp,
         col = "#BBDDEE",
         xlab = "Components",
         main = paste(recountID, featureType, datasetName, "\nvariance barplot"))
    silence <- dev.off()

    ## Plot PC1 vs PC2
    PCplot.file <- file.path(plotDir, paste0(recountID, "_", featureType, "_", datasetName, "_PC1-PC2.pdf"))
    message("\t\tPC 1 vs 2 plot: ", PCplot.file)
    pdf(file = PCplot.file, width = 7, height = 7)
    plot2PCs(dataset, pcs = c(1, 2))
    silence <- dev.off()

    ## Plot PC3 vs PC4
    PCplot.file <- file.path(plotDir, paste0(recountID, "_", featureType, "_", datasetName, "_PC3-PC4.pdf"))
    message("\t\tPC 3 vs 4 plot: ", PCplot.file)
    pdf(file = PCplot.file, width = 7, height = 7)
    plot2PCs(dataset, pcs = c(3, 4))
    silence <- dev.off()

    ## Plot PC5 vs PC6
    PCplot.file <- file.path(plotDir, paste0(recountID, "_", featureType, "_", datasetName, "_PC5-PC6.pdf"))
    message("\t\tPC 5 vs 6 plot: ", PCplot.file)
    pdf(file = PCplot.file, width = 7, height = 7)
    plot2PCs(dataset, pcs = c(5, 6))
    silence <- dev.off()

    ## Combine PC1-PC2 and PC3-PC4  plots in a single figure
    PCplot.file <- file.path(plotDir, paste0(recountID, "_", featureType, "_", datasetName, "_plots.pdf"))
    message("\t\tPC plots: ", PCplot.file)
    pdf(file = PCplot.file, width = 10, height = 10)
    par(mfrow = c(2,2))
    plot(dataset$prcomp, col = "#BBDDEE",
         xlab = "Components",
         main = paste(recountID, featureType, "\n", datasetName, "variance"))
    plot2PCs(dataset, pcs = c(1,2))
    plot2PCs(dataset, pcs = c(3,4))
    plot2PCs(dataset, pcs = c(5,6))
    par(mfrow = c(1,1))
    silence <- dev.off()

  }
}

#### Plot variance per gene at different levels of filtering ####
for (recountID in selectedRecountIDs) {

  parameters <- studyCases[[recountID]]$parameters
  featureType <- parameters$feature
  #  plotFilterHistograms(filteredDataset) #,  plot.file = file.path(parameters$dir$NormalizationImpact, "var_per_gene_hist.pdf"))
  filtered <- studyCases[[recountID]]$datasetsForTest$filtered

  plotFilterHistograms(
    dataset = filtered,
    rawCounts = studyCases[[recountID]]$rawData$countsPerSample,
    #    plot.height = 8,
    plot.file = file.path(
      parameters$dir$NormalizationImpact,
      paste(sep = "_", parameters$recountID, "filtering_variance_per_gene_hist.pdf")))
}

#### Plot library sizes ####
for (recountID in selectedRecountIDs) {
  parameters <- studyCases[[recountID]]$parameters
  featureType <- parameters$feature
  plot.file <- file.path(
    project.parameters$global$dir$figures,
    paste0(recountID, "_", featureType, "_filtered_ranked_libsizes.pdf"))
  message("Plotting ranked library sizes\t", plot.file)
  for (recountID in selectedRecountIDs) {

    parameters <- studyCases[[recountID]]$parameters
    #  plotFilterHistograms(filteredDataset) #,  plot.file = file.path(parameters$dir$NormalizationImpact, "var_per_gene_hist.pdf"))
    filtered <- studyCases[[recountID]]$datasetsForTest$filtered

    ## Add the plot as panel to compare between study cases
    LibsizeRankPlot(count.table = filtered$dataTable,
                    parameters = parameters,
                    plot.file = plot.file)

    # ## Draw a separate plot in the study case-specific directory
    # LibsizeRankPlot(count.table = filtered$dataTable,
    #                 plot.file = file.path(
    #                   parameters$dir$NormalizationImpact,
    #                   paste0(parameters$recountID, "_", featureType, "_filtered_ranked_libsizes.pdf"))
    # )
    # system(paste("open", plot.file))
  }
  # silence <- dev.off()
}

#### Save a memory image that can be re-loaded in the future ####
## This avoids to redo all the normalisation computing each time we want to run the subsequent analyses
featureType <- parameters$feature
if ((project.parameters$global$save.image) &&
    (!project.parameters$global$reload)) {
  featureType <- project.parameters$global$feature
  studyCases.mem.image <- file.path(
    project.parameters$global$dir$memoryImages,
    paste0(
      paste(collapse = "-", selectedRecountIDs),
      "_", featureType,
      "_loaded_studyCase.Rdata"))

  dir.create(project.parameters$global$dir$memoryImages, showWarnings = FALSE, recursive = TRUE)
  message.with.time("Saving memory image after loading: ", studyCases.mem.image)
  save.image(file = studyCases.mem.image)
}


## Indicate that this script has finished running
message.with.time("finished executing 02_load_and_normalise_counts.R")

