#' ###################################################
#' @title Variable ordering according to the p-values returned by differential expression analysis.
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @description Order the variables (columns) of a count table in preparation for supervised
#' classification, by running differential expression analysis with either DESeq2 or edgeR.
#'
#' @param dataTableWithClasses an object containing both the data table and the class label vector.
#' @param method="DESeq2"  choice of method for differential expression analysis. Supported: "DESeq2" , "edgeR".
#' @param edgeRDispEst="tagwise" method used by edgeR to estimate the dispersion. Supported: common, trended, tagwise.
#' @param ... all additional parmeters are passed to the differentail expression method (DEseq2 or edgeR).
#' @return a list with 2 tables:
#' 1) result$geneOrder: a vector of gene names ordered by increasing p-value.
#' 2) result$orderedDataTable ordered count table where columns (genes) have been ordered by increasing p-value.
#' 3) result$DEGtable a table with one row per gene, and one column per DEG statistics (mean, log-ratio, nominal and adjusted p-values, ...)
#'
#' Note that genes are ordered according to nominal p-value (pvalue) rather than
#' adjusted p-value (padj) because DESeq2 produces NA values for the padj.
#'
#' @examples
#'
#' ###################################################
#' ## loading required packages
#' RequiredBioconductorPackages(packages = c("DESeq2","edgeR"))
#'
#' recountID <- "SRP048759"
#' studyCases <- loadCounts(recountID = recountID, mergeRuns = T, classColumn = "tissue")
#' degOrderdPValues <- DEGordering(dataTableWithClasses = studyCases$dataTable, classes =  studyCases$classes, method = "edgeR")
#' ## degPValues<- degPValues[order(degPValues$padj) ,]
#'
#' choose the method of differential expression analysis
#'
#' @export
DEGordering <- function(dataTableWithClasses,
#                       classes,
                       method = "DESeq2",
                       randomized= TRUE,
                       edgeRDispEst="tagwise"){


  RequiredBioconductorPackages(c("DESeq2", "edgeR", "limma"))

  if (!is(dataTableWithClasses, class2 = "DataTableWithClasses")) {
    stop("The function DEGordering() requires an object of class DataTableWithClasses")
  }

  classes <-   dataTableWithClasses$classLabels
  dataTable <-  dataTableWithClasses$dataTable

  result <- list()

  if (method == "DESeq2") {
    message("\t\tInstantiate DESeq2 object ")
    ## Create a DESeqDataset object from the count table
    dds <- DESeqDataSetFromMatrix(dataTable, as.data.frame(classes), ~ classes)

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
    result$geneOrder <- rownames(na.omit(as.data.frame(dataTable)))[geneOrderIndex]
    result$method <- "DESeq2"

    ## Sort the genes (columns) of the count table by increasing p-value
    result$orderedDataTable <- dataTable[result$geneOrder,]
    result$DEGtable <-  DEGtable[result$geneOrder,]
    names(result$DEGtable) <- c(
      "baseMean",
      "log2FC",
      "lfcSE",
      "stat",
      "pvalue",
      "padj" )

    message("\tFinishing from DESeq2 differntial expression analysis")

    if (randomized){
      #result$DEG.DESeq2.randomized <- result$DEG.DESeq2
      result$geneRandomized <-sample( result$geneOrder, replace = F)
      result$randomizedDataTable <- result$orderedDataTable[sample(result$geneRandomized) ,]
      result$randomized <- "randomized-DESeq2"
    }

  } else if (method == "edgeR") {

    ## Build a "model matrix" from the class labels
    message("\t\tInstantiating edgeR object ")

    designMat <- model.matrix(~ classes)
    #dim(designMat)
    #dim(dataTable)

    ## Build dgList object which is required to run edgeR DE analysis
    dgList <- DGEList(counts = (na.omit(as.data.frame(dataTable))))
    #dim(dgList)

    ## Estimate the dispersion parameter for our model. The edgeR method uses
    ## empirical Bayes methods to 'shrink' the genewise dispersion estimates
    ## towards the common dispersion (tagwise dispersion).
    ##
    # Note that either the common or trended dispersion needs to be estimated
    ## before we can estimate the tagwise dispersion.
#    if (edgeRDispEst == "common") {
      message("Estimating common dispersion")
      dgList <- estimateGLMCommonDisp(dgList, design = na.omit(designMat))
#    } else if (edgeRDispEst == "trended") {
      message("Estimating trended dispersion")
      dgList <- estimateGLMTrendedDisp( dgList, design = na.omit(designMat))
#    } else if (edgeRDispEst == "tagwise") {
      message("Estimating tagwise dispersion")
      dgList <- estimateGLMTagwiseDisp(dgList, design = na.omit(designMat))
#    } else {
#      stop(edgedDisp, " is not a valid method for edgeR dispersion estimate. Supported: common, trended, tagwise.")
#    }

    ## Fit edgeR model for differential expression analysis
    message("\t\tedgeR model fitting")
    fit <-   glmFit(dgList, designMat)
#    message("\tFinishig fit settable for closure the fitted object from edgeR ")
    lrt <-  glmLRT(fit)
#   message("\tFinishig lrt settable for closure the fitted object from edgeR ")
    #View(lrt$table)

    ## we can explore the results from topTags function
    ## which give us top DE tags in data frame for a given pair of groups was ranked by p-value or
    ## absolute log-fold change.
    ## edgeRresult<- topTags(lrt)
    ## we ordering the result with the most significance regarding P-Value
    result <- list()
    ## Report the number of NA values for adjusted p-values
    na.padj <- sum(is.na(lrt$table$PValue))
    if (na.padj > 0) {
      message("Beware: edgR reported ", na.padj, " NA for the adjusted p-value.")
      message("we will revome all genes that heve NA values",na.padj,  "to avoid the some problimatics when we passing it for the classifiers")
      #lrt$table$PValue <- na.omit(lrt$table$PValue)
      lrt$table <-na.omit(lrt$table)
    }
    geneOrderIndex <- order(lrt$table$PValue, decreasing = FALSE)
    result$geneOrder <- as.vector(rownames(na.omit(as.data.frame(dataTable)))[geneOrderIndex])

    result$DEGtable <- lrt$table[result$geneOrder,]
    names(result$DEGtable) <- c("log2FC", "logCPM", "LR", "padj")

    result$orderedDataTable <- dataTable[result$geneOrder, ]
    result$method <- "edgeR"
    message("\tFinishing from edgeR differntial expression analysis")


    if (randomized){
      #result$DEG.edgeR.randomized <- result$DEG.edgeR
      result$geneRandomized <- sample( result$geneOrder,replace = F)
      result$randomizedDataTable <- result$orderedDataTable[ sample(result$geneOrder),]
      result$randomized <- "randomaised-edgeR"
    }
    ## end of the if edgeR
  } else {
    stop(method, " is not a valid method for DEGordering. Supported: edgeR, DESeq2")
  }

  return(result)

}
