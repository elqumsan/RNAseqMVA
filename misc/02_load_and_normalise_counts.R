#################### Load counts and pheno ####################
## Load a count Table from recount-experiment, merge counts per sample
## and apply some pre-filtering (remove zero-variance and near-zero-variance genes).

if (project.parameters$global$reload) {
  ## Reload previously stored memory image
  if (!is.null(project.parameters$global$reload.date)) {
    image.date <- project.parameters$global$reload.date
  } else {
    image.date <- Sys.Date()
  }
  studyCases.mem.image <- file.path(
    project.parameters$global$dir$memoryImages,
    paste(sep = "", "loaded_studyCases_",
          paste(collapse = "-", selectedRecountIDs),
          "_", image.date, ".Rdata"))

  message("Reloading study cases from previously stored memory image")
  message("\t", studyCases.mem.image)
  load(studyCases.mem.image)

  ## Reload parameters (their value has been over-written by those loade in memory image)

  ## If requrested, reset the parameters for all the study cases
  ## This is used to re-run the analyses on each study case after
  ## having changed some parameters in the yaml-specific configuration file
  if (project.parameters$global$reload.parameters) {
    project.parameters <- yaml.load_file(configFile)
    project.parameters <- initParallelComputing(project.parameters)
    if (exists("studyCases")) {
      for (recountID in names(studyCases)) {
        parameters <- initRecountID(recountID, project.parameters)
      	studyCases[[recountID]]$parameters <- parameters
	for (dataSetName in names(studyCases[[recountID]]$datasetsForTest)) {
          studyCases[[recountID]]$datasetsForTest[[dataSetName]]$parameters <- parameters
        }
        #  print (studyCases[[recountID]]$parameters$dir$tablesDetail)
      }
    }
  }
  
} else {
  message("Loading study cases")

  studyCases <- list() ## a list containing all the loaded datasets + their pre-processed data
  #recountID <- "SRP057196"
  # recountID <- "SRP056295"
  recountID <- "SRP042620" ## For quick test and debugging
  for (recountID in selectedRecountIDs) {


    message.with.time("Building StudyCase for recountID\t", recountID)

    #### Specify generic and recountID-specific parameters ####
    parameters <- initRecountID(recountID, project.parameters)

    # Main directory should be adapted to the user's configuration
    #  dir.main <- project.parameters$global$dir$main

    # View(parameters)
    studyCases[[recountID]] <- StudyCase(recountID = recountID, parameters = parameters)


    #### Export the count tables with their associated information (pheno table, class labels) in tab-separated value (.tsv) files ###
    exportTables(studyCases[[recountID]])

    ## TO DO LATER: CHECK IF THESE FIGURES ARE WELL GENERATED, AND INCOROPORATE THEM IN THE plot.figure methods

    # plot.file <- file.path(parameters$dir$NormalizationImpact, "log2normCount_hist.pdf")
    # message("\tlog2(norm counts) histogram\t", plot.file)
    # pdf(plot.file, width = 7, height = 5)
    # hist(unlist(studyCases[[recountID]]$log2norm$counts), breaks=100,
    #      col="grey",
    #      main=paste("log2(norm counts) distrib;", recountID),
    #      las=1,
    #      xlab="log2(norm counts)",
    #      ylab="Frequency")
    #
    # silence <- dev.off()

    #   if (ncol(studyCases[[recountID]]$log2norm$dataTable) != length(studyCases[[recountID]]$log2norm$classLabels)){
    #     stop(" the Number of samples in log2norm counts should be the same length of classes")
    #   }
    #
    #
    #
    # } else {
    #   message.with.time("Skipping normalisation for count Table with log2 transformation")
    # }


  } # end loop over recountIDs
  message("Finished loading ", length(studyCases), " study cases")
}


## Compute statistics about loaded datasets
studyCasesStats <- data.frame()
message("Computing summary statistics for ", length(studyCases), " study case(s). ")
for (recountID in selectedRecountIDs) {
  newStats <-
    data.frame(
      recountID = recountID,
      runs = studyCases[[recountID]]$rawData$countsPerRun$nbSamples,
      samples = studyCases[[recountID]]$rawData$countsPerSample$nbSamples,
      samples.filtered = studyCases[[recountID]]$datasetsForTest$filtered$nbSamples,

      genes.ori = studyCases[[recountID]]$rawData$countsPerSample$nbGenes,
      genes.NA = length(studyCases[[recountID]]$datasetsForTest$filtered$naGenes),
      genes.zeroVar = length(studyCases[[recountID]]$datasetsForTest$filtered$zeroVarGenes),
      genes.NZfilter = length(studyCases[[recountID]]$datasetsForTest$filtered$nearZeroVarGenes),
      genes.filtered = studyCases[[recountID]]$datasetsForTest$filtered$nbGenes,

      classes.ori = studyCases[[recountID]]$rawData$countsPerSample$nbClasses,
      classes.filtered = studyCases[[recountID]]$datasetsForTest$filtered$nbClasses
    )

  if (ncol(studyCasesStats) == 0) {
    studyCasesStats <- newStats
  } else {
    studyCasesStats <- rbind(studyCasesStats, newStats)
  }
}
rownames(studyCasesStats) <- studyCasesStats$recountID
studyCasesStats$pc.NA <- 100*studyCasesStats$genes.NA / studyCasesStats$genes.ori
studyCasesStats$pc.zeroVar <- 100*studyCasesStats$genes.zeroVar / studyCasesStats$genes.ori
studyCasesStats$pc.NZfilter <- 100*studyCasesStats$genes.NZfilter / studyCasesStats$genes.ori
studyCasesStats$pc.kept <-  100*studyCasesStats$genes.filtered / studyCasesStats$genes.ori

## TEMPORARY: print out stats about loaded datasets
require(knitr)
kable(t(studyCasesStats))

write.table(x = t(studyCasesStats), file = file.path(parameters$dir$results, "experiment_summaries.tsv"),
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)


gene.pc <- studyCasesStats[, c("pc.NA", "pc.zeroVar", "pc.NZfilter", "pc.kept")]
row.names(gene.pc) <- studyCasesStats$recountID
## Order experiments by increasing percentof kept genes (so that the highest come on top of the barplot)
gene.pc <- gene.pc[order(gene.pc$pc.kept, decreasing = FALSE), ]
## apply(gene.pc, 1, sum) ## The sums must give 100 for each experiment

####  Gene filtering barplots:  proportions of genes filtered out or kept in each study case ####
figure.file <- paste("gene_filtering_summaries.pdf")
barPlot.file <- file.path(project.parameters$global$dir$figures, figure.file)
message("Filtering summary barplot: ", barPlot.file)
pdf(file = barPlot.file, width = 7, height = 2 + 1 * nrow(studyCasesStats))
save.margins <- par("mar")
par(mar = c(5,7,5,1))

kept.label <- paste(sep = "", round(digits = 0, gene.pc$pc.kept), "%")
ypos <- barplot(t(gene.pc), las = 1, horiz = TRUE,
                col = c("black", "red", "orange", "#44DD44"),
                legend.text = c("NA values", "Zero var", "NZ filter", "kept"),
                main = "Filtering impact on study cases",
                xlab = "Proportions of genes",
                xlim = c(0, 170))
text(x = 100, kept.label, y = ypos, pos = 4)

par(mar = save.margins)  ## Restore original graphical parameter
silence <- dev.off()      ## Close graphical device

#### Draw a barplot with the number of samples per class ####
message("exporting barplots with number of samples per class")
for (recountID in names(studyCases)) {
  dir.figures <- file.path(parameters$dir$results, recountID, "figures")
  dir.create(dir.figures, recursive = TRUE, showWarnings = FALSE)
  figure.file <- paste(sep = "", recountID, "_samples_per_classes.pdf")
  barPlot.file <- file.path(dir.figures, figure.file)
  message("\t", recountID, "\tSamples per classes barplot: ", barPlot.file)
  pdf(file = barPlot.file, width = 5, height = 8)
  par.ori <- par(no.readonly = TRUE)
  # par(mfrow = c(4,2))
  classNames <- studyCases[[recountID]]$rawData$countsPerSample$classNames
  classNameLen <- max(nchar(classNames))
  par(mar = c(4, 1 + classNameLen*0.3, 5, 1))
  # par("mar")
  heights <- barplot(
    sort(studyCases[[recountID]]$rawData$countsPerSample$samplesPerClass, decreasing = TRUE),
    horiz = TRUE, las = 1, cex.names = 0.7, main = recountID,
    xlab = "Samples per class", col = "white")
  barplot(sort(studyCases[[recountID]]$datasetsForTest$filtered$samplesPerClass, decreasing = TRUE),
          add = TRUE, horiz = TRUE, las = 1, cex.names = 0.7,
          main = recountID, xlab = "Samples per class",
          col = studyCases[[recountID]]$datasetsForTest$filtered$classColors)
  #          col="#00BB00")
  abline(v = studyCases[[recountID]]$parameters$filtering$minSamplesPerClass, col = "red")
  par(mfrow = c(1,1))
  par(par.ori)
  silence <- dev.off()
}


#### Principal component plots ####
message("\n\tExporting principal component plots")
for (recountID in selectedRecountIDs) {
  message("\tPC plots for study case ", recountID)
  studyCase <- studyCases[[recountID]]
  parameters <- studyCase$parameters

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

    PCplot.file <- file.path(plotDir, paste(sep = "", recountID, "_", datasetName, "_variance.pdf"))
    message("\t\tPC variance plot: ", PCplot.file)
    pdf(file = PCplot.file, width = 7, height = 5)
    plot(dataset$prcomp, col = "#BBDDEE", xlab = "Components", main = paste(recountID, datasetName, "\nvariance barplot"))
    silence <- dev.off()

    ## Plot PC1 vs PC2
    PCplot.file <- file.path(plotDir, paste(sep = "", recountID, "_", datasetName, "_PC1-PC2.pdf"))
    message("\t\tPC 1 vs 2 plot: ", PCplot.file)
    pdf(file = PCplot.file, width = 7, height = 9)
    plot2PCs(dataset, pcs = c(1,2))
    silence <- dev.off()

    ## Plot PC2 vs PC3
    PCplot.file <- file.path(plotDir, paste(sep = "", recountID, "_", datasetName, "_PC3-PC4.pdf"))
    message("\t\tPC 3 vs 4 plot: ", PCplot.file)
    pdf(file = PCplot.file, width = 7, height = 9)
    plot2PCs(dataset, pcs = c(3,4))
    silence <- dev.off()

    ## Combine PC1-PC2 and PC3-PC4  plots in a single figure
    PCplot.file <- file.path(plotDir, paste(sep = "", recountID, "_", datasetName, "_plots.pdf"))
    message("\t\tPC plots: ", PCplot.file)
    pdf(file = PCplot.file, width = 10, height = 10)
    par(mfrow = c(2,2))
    plot(dataset$prcomp, col = "#BBDDEE",
         xlab = "Components", main = paste(recountID, datasetName, "\nvariance barplot"))
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
plot.file <- file.path(
  project.parameters$global$dir$figures,  "filtered_ranked_libsizes.pdf")
message("Plotting ranked library sizes\t", plot.file)
pdf(file = plot.file, width = 7, height = 10)
par(mfrow = c(4, 2))
for (recountID in selectedRecountIDs) {

  parameters <- studyCases[[recountID]]$parameters
  #  plotFilterHistograms(filteredDataset) #,  plot.file = file.path(parameters$dir$NormalizationImpact, "var_per_gene_hist.pdf"))
  filtered <- studyCases[[recountID]]$datasetsForTest$filtered

  ## Add the plot as panel to compare between study cases
  LibsizeRankPlot(count.table = filtered$dataTable, plot.file = NULL)

  ## Draw a separate plot in the study case-specific directory
  LibsizeRankPlot(count.table = filtered$dataTable,
                  plot.file = file.path(
                    parameters$dir$NormalizationImpact,
                    paste(sep = "_", parameters$recountID, "filtered_ranked_libsizes.pdf"))
  )
  # system(paste("open", plot.file))
}
dev.off()

#### Save a memory image that can be re-loaded next time ####
## to avoid re-computing all the normalisation and so on.
if (project.parameters$global$save.image) {
  image.date <- Sys.Date()
  studyCases.mem.image <- file.path(
    project.parameters$global$dir$memoryImages,
    paste(sep = "", "loaded_studyCases_",
          paste(collapse = "-", selectedRecountIDs),
          "_", image.date, ".Rdata"))

  dir.create(project.parameters$global$dir$memoryImages, showWarnings = FALSE, recursive = TRUE)
  message.with.time("Saving memory image after loading: ", studyCases.mem.image)
  save.image(file = studyCases.mem.image)
}


## Indicate that this script has finished running
message.with.time("finished executing 02_load_and_normalise_counts.R")

## Temporary
# source('misc/06_all_variables_vs_all_PCs.R')
