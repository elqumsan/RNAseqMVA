#' @title Apply various filters on the observations and variables of
#' a count table in order to prepare it for classification.
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @description Apply various filters on the observations (e.g. biological samples)
#' and variables (e.g. genes) of
#' a count table in order to prepare it for classification.
#' @param countsWithClasses an object of the class DataTableWithClasses.
#' This object contains the data table + the parameters (including filtering parameters).
#' @param draw.plot=TRUE if TRUE, draw an histogram of variance per gene.
#'
#' @examples
#' ## Load a data set
#' countdata <- loadRecountExperiment(recountID = "SRP048759", mergeRuns = TRUE)
#' dataTable <- t(countdata$merged$sampleCounts)
#' phenoTable <- countdata$merged$samplePheno
#'
#' # Run the filtering
#' filteredData <- filterDataTable(dataTable, phenoTable, classColumn="tissue")
#'
#' # Replace unfiltered data by filtered data
#' dataTable <- filteredData$dataTable
#' phenoTable <- filteredData$phenoTable
#' classLabels <- filteredData$classLabels
#'
#' ## Filter a dataset and build classLabels based on 2 columns
#' countdata <- loadRecountExperiment(recountID = "SRP057196", mergeRuns = TRUE)
#' dataTable <- t(countdata$merged$sampleCounts)
#' phenoTable <- countdata$merged$samplePheno
#'
#' #' # Run the filtering
#' filteredData <- filterDataTable(dataTable, phenoTable, classColumn=c("tissue", "cell.type"), minSamplesPerClass=5)
#' table(filteredData$classLabels)
#'
#' @import caret
#' @import lattice
#' @import ggplot2
#' @export
filterDataTable <- function(rawCounts,
                            draw.plot = TRUE) {

  RequiredBioconductorPackages(c("recount", "SummarizedExperiment"))


  message.with.time("Filtering count table")

  ## Check validity of the input class
  if (!is(rawCounts, "DataTableWithClasses")) {
    stop("filterDataTable() requires an object of class DataTableWithClasses")
  }


  ## Take parameters from the DataTableWithClasses object
  parameters <- rawCounts$parameters

  ## Check required parameters
  if (is.null(parameters$filtering$na.rm)) {
    message("\tfilterDataTable()\tundefined filtering parameter: na.rm. Setting to default (TRUE).")
    parameters$filtering$na.rm <- TRUE
  }
  na.rm <- parameters$filtering$na.rm

  if (is.null(parameters$filtering$na.rm)) {
    message("\tfilterDataTable()\tundefined filtering parameter: minSamplesPerClass. Setting to default (10).")
    parameters$filtering$minSamplesPerClass <- 10
  }
  minSamplesPerClass <- parameters$filtering$minSamplesPerClass

  if (is.null(parameters$filtering$nearZeroVarFilter)) {
    message("\tfilterDataTable()\tundefined filtering parameter: nearZeroVarFilter. Setting to default (FALSE).")
    parameters$filtering$nearZeroVarFilter <- FALSE
  }
  nearZeroVarFilter <- parameters$filtering$nearZeroVarFilter

  ## Set the updated parameters to the object
  rawCounts$parameters <- parameters

  # for (p in c("na.rm", "minSamplesPerClass", "nearZeroVarFilter")) {
  #   if (is.null(parameters[[p]])) {
  #     stop("Missing required parameter: '", p,
  #          "'.\n\tPlease check configuration file. ")
  #   } else {
  #     assign(p, parameters[[p]])
  #   }
  # }



  ## Initialise a vector with the genes that will be kept through the different filtering stpes.
  keptGenes <- rawCounts$geneNames


  ## Check if there are NA values, and discard all genes having at least one NA value
  if (sum(is.na(rawCounts$dataTable)) > 0) {
    naGenes <- rownames(rawCounts$dataTable)[apply(is.na(rawCounts$dataTable), 1, sum) > 0]
    ## Report the number of genes with NA values
    message("\tThis count Table contains ", length(naGenes), " genes with NA values")
    message("\tFiltering out ", length(naGenes)," genes with NA values")
  } else {
    message("\tThis count Table does not contain any NA values. ")
    naGenes <- vector()
  }

  if (na.rm) {
    keptGenes <- setdiff(rawCounts$geneNames, naGenes)
    message("\tNA filter: discarding ", length(naGenes), " genes out of ", rawCounts$nbGenes, "; keeping ", length(keptGenes))
  }
  #filteredDataTable <- rawCounts$dataTable[, keptGenes]


  ## Detect genes with zero variance
  message("\tDetecting genes with zero variance")
  varPerGene <- apply(rawCounts$dataTable, 1 , var, na.rm=TRUE) # Compute variance per gene (row)
  # table(varPerGene == 0)

  # sum(zeroVarGenes == 0) # count genes with zero variance
  zeroVarGenes <- names(varPerGene)[varPerGene == 0] # identify genes having zero variance
  # length(zeroVarGenes) ## Number of genes to discard
  keptGenes <- setdiff(keptGenes, zeroVarGenes)
  #length(keptGenes) ## Number of genes to keep
#  filteredDataTable <- rawCounts$dataTable[keptGenes, ]
  message("\tZero var filter: discarding ", length(zeroVarGenes)," genes with zero variance; keeping ", length(keptGenes), " genes")
  # dim(countsWithClasses$dataTable)
  # dim(filteredDataTable)



  ## Use caret::nearZeroVar() to discard supposedly bad predictor genes.
  ## This includes zero variance genes (already filtered above) but also less trivial cases, see nearZeroVar() doc for details.
  nearZeroVarGenes <- NULL
  if (nearZeroVarFilter) {
    message("\tDetecting genes with near-zero variance using caret::nearZeroVar()")
    nearZeroVarIndices <- nearZeroVar(t(rawCounts$dataTable), allowParallel = TRUE, saveMetrics = FALSE)
    nearZeroVarGenes <- rownames(rawCounts$dataTable)[nearZeroVarIndices]
    nearZeroVarGenes <- setdiff(nearZeroVarGenes, zeroVarGenes)
    keptGenes <- setdiff(keptGenes, nearZeroVarGenes)
    message("\tNear zero var filter: discarding ", length(nearZeroVarGenes)," genes with near-zero variance (poor predictors); kept genes: ", length(keptGenes))
  }

  ## Extract the table of counts with the subset of kept genes at the end of the different gene filters.
  filteredDataTable <- rawCounts$dataTable[keptGenes,]
  # dim(filteredDataTable)

  ## Plot an histogram to compare variance distribution  between all genes and those with  near-zero variance
  if (draw.plot) {
    plot.file <- file.path(parameters$dir$NormalizationImpact, "var_per_gene_hist.pdf")
    message("\tVariance per gene histograms\t", plot.file)
    pdf(plot.file, width=7, height=12)

    logVarPerGene <- log2(varPerGene)
    noInfVar <- !is.infinite(logVarPerGene)
    xmin <- floor(min(logVarPerGene[noInfVar]))
    xmax <- ceiling(max(logVarPerGene[noInfVar]))
    xlim <- c(xmax, xmin)
    varbreaks <- seq(from=xmin,  to=xmax, by=0.1)
    if (nearZeroVarFilter) {
      par(mfrow=c(4,1))
    } else {
      par(mfrow=c(3,1))
    }
    hist(log2(varPerGene[noInfVar]),
         breaks=varbreaks,
         col="gray", border = "gray",
         main = paste("All non-zero var genes; ", parameters$recountID),
         xlab="log2(varPerGene)",
         ylab="Number of genes")
#    legend("topright", parameters$recountID)


    if (nearZeroVarFilter) {
      hist(log2(varPerGene[nearZeroVarGenes]),
         breaks=varbreaks,
         col="red", border = "orange",
         main = "Near zero variance",
         xlab="log2(varPerGene)",
         ylab="Number of genes")
    }
    hist(log2(varPerGene[keptGenes]),
         breaks=varbreaks,
         col="darkgreen", border = "#00BB00",
         main = "Kept genes",
         xlab="log2(varPerGene)",
         ylab="Number of genes")

    ## Count the number of zero values per gene
    zerosPerGene <- apply(rawCounts$dataTable == 0, 1, sum)
    zerobreaks <- seq(from=0, to=max(zerosPerGene+1), length.out =50)
    # zerobreaks <- seq(from=0, to=max(zerosPerGene+1), by=1)

    #### Histogram of zero values per gene. ####
    ##
    ## Displayed in green (color for kept genes) because we will then
    ## overlay the histograms of near zero var (orange) and zero var (red),
    ## so that the remaining part of the histogram will correspond to genes
    ## kept.
    hist(zerosPerGene,
         breaks = zerobreaks,
         main = "Zeros per gene",
         xlab = "Number of zero values",
         ylab = "Number of genes",
         col = "#00BB00", border = "#00BB00")
    if (nearZeroVarFilter) {
      hist(zerosPerGene[union(nearZeroVarGenes, zeroVarGenes)],
         breaks = zerobreaks,
         add = TRUE, col="orange", border="orange")
    }
    hist(zerosPerGene[zeroVarGenes],
         breaks = zerobreaks,
         add = TRUE, col = "red", border = "red")
    if (nearZeroVarFilter) {
      legend("topleft",
             legend=paste(
               sep= ": ",
               c("Kept genes", "Near-zero variance", "Zero variance"),
               c(length(keptGenes), length(nearZeroVarGenes), length(zeroVarGenes))),
             lwd=5,
             col=c("#00BB00", "orange", "red")
      )
    } else {
      legend("topleft",
             legend=paste(
               sep= ": ",
               c("Kept genes",  "Zero variance"),
               c(length(keptGenes), length(zeroVarGenes))),
             lwd=5,
             col=c("#00BB00", "red"))

    }
    #    hist(zerosPerGene[nearZeroVarGenes])

    par(mfrow=c(1,1))

    silence <- dev.off()

  }


  ################################################################
  ## Check if there are NA values in the sample classes
  discardedSamples <- is.na(rawCounts$classLabels)
  if (sum(discardedSamples) > 0) {
    message("\tDiscarding ", sum(discardedSamples), " samples with undefined class in ", classColumn, " column of pheno table")
    nonaSamples <- !is.na(rawCounts$classLabels)
    filteredDataTable <- filteredDataTable[, nonaSamples]
    # dim(filteredDataTable)
    # dim(filteredPhenoTable)
    filteredPhenoTable<- rawCounts$phenoTable[nonaSamples, ]
    filteredClassLabels <- rawCounts$classLabels[nonaSamples]
    # table(countsWithClasses$classLabels, filteredClassLabels)
  } else {
    filteredPhenoTable<- rawCounts$phenoTable
    filteredClassLabels <- rawCounts$classLabels
  }

  ################################################################
  ## Select classess for which we dispose of at least 10 samples
  # length(countsWithClasses$classLabels)
  message("\tSelecting classes with at least ", minSamplesPerClass, " samples")
  samplesPerClass <- table(filteredClassLabels)
  discardedClasses <- names(samplesPerClass)[samplesPerClass < minSamplesPerClass]
  keptClasses <- names(samplesPerClass)[samplesPerClass >= minSamplesPerClass]
  keptSamples <- filteredClassLabels %in% keptClasses
  if (length (discardedClasses) > 0) {
    message("\tDiscarding ", length (discardedClasses), " classes containing less than ", minSamplesPerClass, " samples")
    message("\tDiscarded classes\t", paste(collapse=",", discardedClasses))
  }
  message("\tKeeping ", sum(keptSamples), " samples from ",
          length(keptClasses), " classes")
  message("\tKept classes\t", paste(keptClasses, collapse=", "))


  ## Update count table, pheno table and countsWithClasses$classLabels vector to keep only the samples belonging to selected classess
  filteredDataTable <- filteredDataTable[, keptSamples]
  # dim(filteredDataTable)
  filteredPhenoTable <- filteredPhenoTable[keptSamples, ]
  # dim(filteredPhenoTable)
  filteredClassLabels <- filteredClassLabels[keptSamples]
  # length(filteredClassLabels)
  filteredSampleNb <- nrow(filteredDataTable)


  ## Return unfiltered count table + phenotable + countsWithClasses$classLabels
  result <- DataTableWithClasses(dataTable= filteredDataTable,
                                  phenoTable =  filteredPhenoTable,
                                  # classColumn = rawCounts$classColumn,
                                  # classColors = rawCounts$classColors,
                                  # variablesType = "all",
                                  dataType = "filtered_counts",
                                  parameters = rawCounts$parameters)

  # class(result)
  summary(result)

  # result$dataTable <- filteredDataTable
  # result$phenoTable <- filteredPhenoTable
  # result$classLabels <- filteredClassLabels
  # result$nbSamples <- ncol(filteredDataTable)
  # result$nbGenes <- nrow(filteredDataTable)

  ## Include the filtering criteria in the returned list
  result$naGenes <- naGenes
  result$zeroVarGenes <- zeroVarGenes
  result$nearZeroVarGenes <- nearZeroVarGenes
  result$keptGenes <- keptGenes
  result$keptClasses <- keptClasses

  message.with.time("Finished filterDataTable() for Recount experiment ID ", parameters$recountID)
  ################################################################################
  ## Return an object with the filtered counts + updated pheno table selected classes and samples
  ## Transpose the table in order to get it in the suitable format for classifiers:
  ## one row per individual, one column per variable.
  message("\tBefore filtering: ",
          rawCounts$nbGenes, " genes x ",
          rawCounts$nbSamples, " samples",
          " belonging to ", rawCounts$nbClasses, " classes")
  message("\tAfter filtering: ",
          result$nbGenes, " genes x ",
          result$nbSamples, " samples",
          " belonging to ", result$nbClasses, " classes")
  # View(filteredDataTable)
  # dim(filteredDataTable)


  return(result)

}
