
#' @title plot histograms of the variance distributions at different steps of the filtering process
#' @author Jacques van Helden
#' @param dataset an object of class DatasetWithClasses returned by filterDataTable()
#' @param plot.file=NULL save the plot in a specified pdf file
#' @export
plotFilterHistograms <- function(dataset,
                                 rawCounts,
                                 plot.file = NULL,
                                 plot.height=NULL) {

  message("\t", dataset$ID, "\tHistograms of variance per feature\t", plot.file)

  ## Check the class of input object
  if (!is(dataset, "DataTableWithClasses")) {
    stop("plotFilterHistograms()\tdataset should belong to class DataTableWithClasses. ")
  }

  ## Get variance per gene
  if (is.null(dataset$varPerGeneRaw)) {
    stop("plotFilterHistograms()\tdataset must contain a field varPerGeneRaw as computed by filterDataTable()")
  }
  varPerGene <- na.omit(dataset$varPerGeneRaw)
  # summary(varPerGene)

  ## Get list of kept genes
  if (is.null(dataset$keptGenes)) {
    stop("plotFilterHistograms()\tdataset must contain a field keptGenes as computed by filterDataTable()")
  }
  keptGenes <- dataset$keptGenes

  ## Get list of genes with zero variance
  if (is.null(dataset$zeroVarGenes)) {
    stop("plotFilterHistograms()\tdataset must contain a field zeroVarGenes as computed by filterDataTable()")
  }
  zeroVarGenes <- dataset$zeroVarGenes

  ## Near-zero filter
  parameters <- dataset$parameters
  featureType <- parameters$feature
  if (is.null(parameters$filtering$nearZeroVarFilter)) {
    nearZeroVarFilter <- FALSE
    nb.panels <- 3
  } else {
    nearZeroVarFilter <- parameters$filtering$nearZeroVarFilter
    nb.panels <- 4
  }

  ## Open pdf file if required
  if (!is.null(plot.file)) {
    if (is.null(plot.height)) {
      plot.height <- nb.panels * 2.5
    }
    pdf(plot.file, width = 7, height = plot.height)
  }


  logVarPerGene <- log2(varPerGene)
  noInfVar <- !is.infinite(logVarPerGene)
  xmin <- floor(min(logVarPerGene[noInfVar]))
  xmax <- ceiling(max(logVarPerGene[noInfVar]))
  xlim <- c(xmin, xmax)
  varbreaks <- seq(from = xmin,  to = xmax, by  = 0.1)
  if (nearZeroVarFilter) {
    par(mfrow = c(4,1))
  } else {
    par(mfrow = c(3,1))
  }
  hist(log2(varPerGene[noInfVar]),
       breaks = varbreaks,
       col = "gray", border = "gray",
       main = paste("All non-zero var genes;", parameters$recountID, featureType, parameters$short_label),
       xlab = "log2(varPerGene)",
       ylab = "Number of genes",
       xlim = xlim)
  #    legend("topright", parameters$recountID)

  ## Plot variance histogram for genes filtered out by the near-zero variance filter
  if (nearZeroVarFilter) {
    ## Get list of near-zero variance genes
    if (is.null(dataset$nearZeroVarGenes)) {
      stop("plotFilterHistograms()\tdataset must contain a field nearZeroVarGenes as computed by filterDataTable()")
    }
    nearZeroVarGenes <- dataset$nearZeroVarGenes


    hist(log2(varPerGene[nearZeroVarGenes]),
         breaks = varbreaks,
         col = "red", border = "orange",
         main = paste("Near zero variance; ", parameters$recountID, featureType, parameters$short_label),
         xlab = "log2(varPerGene)",
         ylab = "Number of genes",
         xlim = xlim)
  }

  ## Plot variance histogram for genes kept after filtering
  hist(log2(varPerGene[keptGenes]),
       breaks = varbreaks,
       col = "darkgreen", border = "#00BB00",
       main = paste("Kept genes; ", parameters$recountID, featureType, parameters$short_label),
       xlab = "log2(varPerGene)",
       ylab = "Number of genes",
       xlim = xlim)

  ## Count the number of zero values per gene
  zerosPerGene <- na.omit(apply(rawCounts$dataTable == 0, 1, sum))
  # summary(zerosPerGene)
  zerobreaks <- seq(from = 0, to = max(zerosPerGene + 1))
  # zerobreaks <- seq(from=0, to=max(zerosPerGene+1), by=1)

  #### Histogram of zero values per gene. ####
  ##
  ## Displayed in green (color for kept genes) because we will then
  ## overlay the histograms of near zero var (orange) and zero var (red),
  ## so that the remaining part of the histogram will correspond to genes
  ## kept.
  hist(zerosPerGene,
       breaks = zerobreaks,
       main = paste("Zeros per gene; ", parameters$recountID, featureType, parameters$short_label),
       xlab = "Number of zero values",
       ylab = "Number of genes",
       col = "#00BB00", border = "#00BB00")
  if (nearZeroVarFilter) {
    hist(zerosPerGene[union(nearZeroVarGenes, zeroVarGenes)],
         breaks = zerobreaks,
         add = TRUE, col = "orange", border = "orange")
  }
  hist(zerosPerGene[zeroVarGenes],
       breaks = zerobreaks,
       add = TRUE, col = "red", border = "red")
  if (nearZeroVarFilter) {
    legend("top",
           legend = paste(
             sep = ": ",
             c("Kept genes", "Near-zero variance", "Zero variance"),
             c(length(keptGenes), length(nearZeroVarGenes), length(zeroVarGenes))),
           lwd = 5,
           cex = 0.8,
           col = c("#00BB00", "orange", "red")
    )
  } else {
    legend("top",
           legend = paste(
             sep = ": ",
             c("Kept genes",  "Zero variance"),
             c(length(keptGenes), length(zeroVarGenes))),
           lwd = 5,
           cex = 0.8,
           col = c("#00BB00", "red"))

  }
  #    hist(zerosPerGene[nearZeroVarGenes])

  par(mfrow = c(1,1))
  if (!is.null(plot.file)) {
    silence <- dev.off(); rm(silence)
  }
}
