#################### Load counts and pheno ####################
## Load a count Table from recount-experiment, merge counts per sample
## and apply some pre-filtering (remove zero-variance and near-zero-variance genes).
studyCases.mem.image <- file.path(
  paste(project.parameters$global$dir$memoryImages,
  "in_DATE",Sys.Date(),"loaded_studyCases.Rdata",sep = "_"))


if (project.parameters$global$reload) {
  ## Reload previously stored memory image
  message("Reloading study cases from previously stored memory image")
  message("\t", studyCases.mem.image)
  load(studyCases.mem.image)
} else {
  message("Loading study cases")

  studyCases <- list() ## a list containing all the loaded datasets + their pre-processed data

  #recountID <- "SRP056295"
  recountID <- "SRP042620" ## For quick test and debugging
  for (recountID in selectedRecountIDs) {


    message.with.time ("Building StudyCase for recountID\t", recountID)

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
    # pdf(plot.file, width=7, height=5)
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

figure.file <- paste("experiments_summaries.pdf")
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

par(mar = save.margins)  ## Restore original graphical parameter
silence<- dev.off()      ## Close graphical device

#### Draw a barplot with the number of samples per class ####
figure.file <- paste("samples_per_classes.pdf")
barPlot.file <- file.path(parameters$dir$results,figure.file)
message("Samples per classes barplot: ", barPlot.file)
pdf(file = barPlot.file, width=8, height=12)
save.margins <- par("mar")
par(mfrow=c(4,2))
par(mar=c(4,15,5,1))
for (recountID in names(studyCases)) {
  heights <- barplot(
    sort(studyCases[[recountID]]$rawData$countsPerSample$samplesPerClass, decreasing = TRUE),
    horiz = TRUE, las=1, cex.names = 0.7, main=recountID,
    xlab="Samples per class", col="white")
  barplot(sort(studyCases[[recountID]]$datasetsForTest$filtered$samplesPerClass, decreasing = TRUE),
          add=TRUE, horiz = TRUE, las=1, cex.names = 0.7,
          main=recountID, xlab="Samples per class",
          col = studyCases[[recountID]]$datasetsForTest$filtered$classColors)
#          col="#00BB00")
  abline(v = studyCases[[recountID]]$parameters$minSamplesPerClass, col="red")
}
par(mfrow=c(1,1))
par(mar = save.margins)
silence<- dev.off()

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


#### Principal component plots ####
for (recountID in selectedRecountIDs) {
  plotDir <- studyCases[[recountID]]$parameters$dir$PCviz
  message("plotDir = ", plotDir)

  ## Variance barplot of the components
  PCplot.file <- file.path(plotDir, paste(sep="", recountID, "_log2norm_prcomp_variance.pdf"))
  message("PC variance plot: ", PCplot.file)
  pdf(file = PCplot.file, width=7, height=5)
  plot(studyCases[[recountID]]$datasetsForTest$log2normPCs$prcomp, col="#BBDDEE", xlab="Components")
  silence <- dev.off()

  #biplot(studyCases[[recountID]]$datasetsForTest$log2normPCs$prcomp, pc.biplot = TRUE)

  ## Plot PC1 vs PC2
  PCplot.file <- file.path(plotDir, paste(sep="", recountID, "_log2norm_PC1-PC2.pdf"))
  message("PC plot: ", PCplot.file)
  pdf(file = PCplot.file, width=7, height=9)
  plot2PCs(studyCases[[recountID]]$datasetsForTest$log2normPCs, pcs = c(1,2))
  silence <- dev.off()

  ## Plot PC2 vs PC3
  PCplot.file <- file.path(plotDir, paste(sep="", recountID, "_log2norm_PC3-PC4.pdf"))
  message("PC plot: ", PCplot.file)
  pdf(file = PCplot.file, width=7, height=9)
  plot2PCs(studyCases[[recountID]]$datasetsForTest$log2normPCs, pcs = c(3,4))
  silence <- dev.off()

  ## Combine PC1-PC2 and PC3-PC4  plots in a single figure
  PCplot.file <- file.path(plotDir, paste(sep="", recountID, "_log2norm_PCplots.pdf"))
  message("PC plot: ", PCplot.file)
  pdf(file = PCplot.file, width=10, height=10)
  par(mfrow=c(2,2))
  plot(studyCases[[recountID]]$datasetsForTest$log2normPCs$prcomp, col="#BBDDEE",
       xlab="Components", main=paste(recountID, " PC variance baplot"))
  plot2PCs(studyCases[[recountID]]$datasetsForTest$log2normPCs, pcs = c(1,2))
  plot2PCs(studyCases[[recountID]]$datasetsForTest$log2normPCs, pcs = c(3,4))
  plot2PCs(studyCases[[recountID]]$datasetsForTest$log2normPCs, pcs = c(5,6))
  par(mfrow=c(1,1))
  silence <- dev.off()
}


## Save a memory image that can be re-loaded next time
## to avoid re-computing all the normalisation and so on.
if (project.parameters$global$save.image) {
  dir.create(project.parameters$global$dir$memoryImages, showWarnings = FALSE, recursive = TRUE)
  message.with.time("Saving memory image after loading: ", studyCases.mem.image)
  save.image(file = studyCases.mem.image)
}

## Indicate that this script has finished running
message.with.time("finished executing 02_load_and_normalise_counts.R")


