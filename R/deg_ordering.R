#' ###################################################
#' @title Variable ordering according to the p-values returned by differential expression analysis.
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @description Order the variables (columns) of a count table in preparation for supervised
#' classification, by running differential expression analysis with either DESeq2 or edgeR.
#'
#' @param countDataset an object of class DataTableWithClasses containing counts of reads per features. Importantly, edgeR and DESeq2 require non-normalized counts as input.
#' @param permuteLabels=FALSE if TRUE, class labels are permuted before running the differential analysis, in order to run a negative control (empirical approximation of the null hypothesis)
#' @param method="DESeq2"  choice of method for differential expression analysis. Supported: "DESeq2" , "edgeR".
#' @param alpha=0.05 threshold on adjusted p-value (FDR) to call genes positive
#' @param ... all additional parmeters are passed to the differential expression method (DEseq2 or edgeR).
#'
#' @return a list with the following fields.
#' \itemize{
#' \item method: method specified in the function call
#' \item permuteLabels: Boolean variable indicating whether the test was led with permuted class labels (negative control)
#' \item classLabels: vector with the class labels used for the analysis (the ones from the original data, or the permuted ones if the option permuteLabels was TRUE)
#' \item geneOrder: a vector of gene names ordered by increasing p-value.
#' \item DEGtable a table with one row per gene, and one column per DEG statistics (mean, log-ratio, nominal and adjusted p-values, ...)
#' \item orderedDataTable ordered count table where rows (genes) have been ordered by increasing p-value.
#' }
#' Note that genes are ordered according to nominal p-value (pvalue) rather than
#' adjusted p-value (padj) because DESeq2 produces NA values for the padj.
#'
#' @examples
#'
#' ###################################################
#' ## loading required packages
#'
#' recountID <- "SRP048759"
#' studyCases <- loadCounts(recountID = recountID, mergeRuns = T, classColumn = "tissue")
#' filteredCounts <- studyCases[[recountID]]$datasetsForTest$filtered
#' degOrderdPValues <- DEGordering(countDataset = filteredCounts, method = "edgeR")
#' ## degPValues<- degPValues[order(degPValues$padj) ,]
#'
#' @export
DEGordering <- function(countDataset,
                        method = "DESeq2",
                        permuteLabels = FALSE,
                        alpha = 0.05) {

  LoadRequiredBioconductorPackages(c("DESeq2", "edgeR"))

  #### Check parameters ####
  if (!is(countDataset, class2 = "DataTableWithClasses")) {
    stop("The function DEGordering() requires an object of class DataTableWithClasses")
  }

  ## Instantiate the result object
  result <- list()

  ## Include the differential analysis method in the result
  result$method <- method

  ## Define class labels (permuted or not) and include them in the result
  classLabels <-   countDataset$classLabels
  if (permuteLabels) {
    classLabels <- sample(x = classLabels, replace = FALSE)
  }
  result$classLabels <- classLabels ## include class labels in the result
  # table(classLabels)
  # length(unique(classLabels))

  message("\tDEGordering()")
  message("\t\t", nrow(countDataset$dataTable), " genes\t")
  message("\t\t", ncol(countDataset$dataTable), " samples\t")
  message("\t\t", length(unique(classLabels)), " classes\t")


  if (method == "DESeq2") {

    #### DESeq2 analysis ####
    message("\t\tInstantiate DESeq2 object ")
    ## Create a DESeqDataset object from the count table
    dds <- DESeqDataSetFromMatrix(countDataset$dataTable, as.data.frame(classLabels), ~ classLabels)

    ## Run  differential expression analysis with DESeq2
    dds <- DESeq(dds)

    # Cast the results from DESeq differential expression analysis
    DEGtable <- results(dds, independentFiltering = F )

    ## Report the number of NA values for adjusted p-values
    na.padj <- sum(is.na(DEGtable$padj))
    if (na.padj > 0) {
      message("Beware: DESeq2 reported ", na.padj, " NA for the adjusted p-value.")
      message("we will revome all genes that heve NA values",na.padj,  "to avoid the some problimatics when we passing it for the classifiers")
      #DEGtable$padj <- na.omit( DEGtable$padj)
      DEGtable <- na.omit( DEGtable)

      # filtered.DEGtable <- as.vector(na.omit( as.data.frame( DEGtable$padj)))
      # truepadj <- DEGtable$padj[filtered.DEGtable]
      # rownames(filtered.DEGtable) <-rownames(DEGtable)[DEGtable$padj == filtered.DEGtable ]
    }
    #as.data.frame(DEGtable)
    ## Notes: we have some NA values for adjusted p-value so that we will choose the p Value for ordering
    geneOrderIndex <- order(na.omit( DEGtable$padj), decreasing = FALSE)
    result$geneOrder <- rownames(na.omit(as.data.frame(countDataset$dataTable)))[geneOrderIndex]
    #    result$method <- "DESeq2"

    ## Sort the genes (columns) of the count table by increasing p-value
    result$orderedDataTable <- countDataset$dataTable[result$geneOrder,]
    result$DEGtable <-  DEGtable[result$geneOrder,]
    names(result$DEGtable) <- c(
      "baseMean",
      "log2FC",
      "lfcSE",
      "stat",
      "pvalue",
      "padj" )

    message("\tCompleted differential expression analysis with DESeq2")

    # if (randomized){
    #   #result$DEG.DESeq2.randomized <- result$DEG.DESeq2
    #   result$geneRandomized <-sample(result$geneOrder, replace = FALSE)
    #   result$randomizedDataTable <- result$orderedDataTable[sample(result$geneRandomized) ,]
    #   result$randomized <- "randomized-DESeq2"
    # }


  } else if (method == "edgeR") {
    #### edgeR analysis ####
    message("\t\tInstantiating edgeR object ")

    ## Build a "model matrix" from the class labels
    ## This matrix contains one row per sample and one column per class
    designMat <- model.matrix(~ classLabels)
    # dim(designMat)
    # dim(countDataset$dataTable)
    # View(designMat)

    ## Build dgList object which is required to run edgeR DE analysis
    dgList <- DGEList(counts = countDataset$dataTable)
    # names(dgList)
    # dim(dgList)
    # class(dgList)
    # View(dgList)

    ## Estimate the dispersion parameter for our model. The edgeR method uses
    ## empirical Bayes methods to 'shrink' the genewise dispersion estimates
    ## towards the common dispersion (tagwise dispersion).
    ##
    # Note that either the common or trended dispersion needs to be estimated
    ## before we can estimate the tagwise dispersion.
    message("\t\tEstimating dispersion")
    dgList <- estimateDisp(dgList, design = designMat)

    # ## Fit edgeR model for differential expression analysis.
    ## We chose glmQLFit because it is claimed to offer a more accurate control of type I error
    message("\t\tedgeR model fitting with glmQLFit()")
    fit <- glmQLFit(dgList, designMat)
    # #    message("\tFinishig fit settable for closure the fitted object from edgeR ")
    # lrt <-  glmLRT(fit)
    # #   message("\tFinishig lrt settable for closure the fitted object from edgeR ")
    # View(lrt$table)

    ## Run test to detect differentially expressed genes
    qlf <- glmQLFTest(fit, coef=2:ncol(designMat))
    # class(qlf)
    # names(qlf)
    qlf.TT <- topTags(qlf, n = nrow(qlf$table))
    # dim(qlf.TT)
    # names(qlf.TT)
    # View(as.data.frame(qlf.TT))

    ## Select genes called positive

    ## Plot multiple MA plots, since we have one logFC for each condition relative to the first one
    # i <- 2 # for quick test
    for (i in 2:ncol(designMat)) {
      logFC.column <- i - 1
      condition <- names(qlf.TT$table)[logFC.column]
      logFC <- qlf.TT$table[,logFC.column]
      logCPM <- qlf.TT$table$logCPM
      plot(logCPM, logFC,
           pch = 1,
           xlab = "logCPM",
           ylab = "logFC",
           col = densCols(x = logCPM, y = logFC),
           main = paste(sep = "", "condition ", i, " vs 1"))
      grid()
      abline(h = 0, col = "black")
    }
    # topTags(qlf)
    # names(topTags(qlf))
    # names(qlf$table)
    # head(qlf$table)
    # View(qlf$table)

    with(topTags(qlf), plot(logCPM, PValue, pch=16, cex=0.2, log="y"))

    ## we can explore the results from topTags function
    ## which give us top DE tags in data frame for a given pair of groups was ranked by p-value or
    ## absolute log-fold change.
    ## edgeRresult<- topTags(lrt)
    ## we ordering the result with the most significance regarding P-Value
    #    result <- list()
    ## Report the number of NA values for adjusted p-values
    na.padj <- sum(is.na(lrt$table$PValue))
    if (na.padj > 0) {
      message("Beware: edgR reported ", na.padj, " NA for the adjusted p-value.")
      message("we will revome all genes that heve NA values",na.padj,  "to avoid the some problimatics when we passing it for the classifiers")
      #lrt$table$PValue <- na.omit(lrt$table$PValue)
      lrt$table <-na.omit(lrt$table)
    }
    geneOrderIndex <- order(lrt$table$PValue, decreasing = FALSE)
    result$geneOrder <- as.vector(rownames(na.omit(as.data.frame(countDataset$dataTable)))[geneOrderIndex])

    result$DEGtable <- lrt$table[result$geneOrder,]
    names(result$DEGtable) <- c("log2FC", "logCPM", "LR", "padj")

    result$orderedDataTable <- countDataset$dataTable[result$geneOrder, ]
    #result$method <- "edgeR"
    message("\tFinished edgeR differntial expression analysis")


    # if (randomized){
    #   #result$DEG.edgeR.randomized <- result$DEG.edgeR
    #   result$geneRandomized <- sample( result$geneOrder,replace = F)
    #   result$randomizedDataTable <- result$orderedDataTable[ sample(result$geneOrder),]
    #   result$randomized <- "randomaised-edgeR"
    # }
    ## end of the if edgeR
  } else {
    stop(method, " is not a valid method for DEGordering. Supported: edgeR, DESeq2")
  }

  return(result)

}
