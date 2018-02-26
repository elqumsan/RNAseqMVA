#' @title study the effect of normalisation onto the efficiency of the classifier
#' @author Mustafa ABUELQUMSAN and Jacques van Helden
#' @description this script measures the impact of normalisation to improve the effeciancy
#' of classification. We start by calculateing the miscalssification errors when the classifier
#' is applied to the raw counts directly without any prior treatment. We then perform some
#' "Normalisation" treatment onto these raw count and study the miscalssification errors
#' in order to evaluate if normalization improves the accuracy.
#'
#' @param countTable
#' @param classes
#' @param classfier
#' @param trainingproportion
#' @param k # number of neighbour for the classfier


## TO DO (next steps)
##
## - study the impact of parameterson classifiers. IN particular, k for KNN
## - separate the question about variable number and variable ordering -> analyse data with increasing number of variables (genes) picked up at random, but with the real class labels.
## - add SVM to the list of classifiers
## - clustering and heat map
## - DEG graphics: Volcano plot, MA plot, p-value distribution
## - DEG: negative control with permuted classes


## dir.knnSubFigures<-file.path(dir.figures, paste("Type_of_classifier_knn" ,sep=""))
## dir.create(dir.knnSubFigures, showWarnings = F, recursive = T)
## dir.rfSubFigures <-file.path(dir.figures, paste("Type_of_classifier_rf",sep = ""))
## dir.create(dir.rfSubFigures, showWarnings = F, recursive = T)

#############################################################################################
###### Draw some plots to showing effects of some treatments(Norm and log2) on the Count Table
#############################################################################################
message.with.time("Drawing plots describing count table statistics")

## histogram of normalized counts. Zoom on the representative part of the histogram
#setEPS()
#postscript(paste(parameters$recountID, "counts_norm-perc75_hist", ".eps", sep = ""))
#postscript(paste(" The normalisation count", ".eps", sep = ""))
pdf(file=file.path(
  dir.NormImpact,
  paste(parameters$recountID, "_counts_norm-perc75_hist", ".pdf", sep = "")),
  width = 8, height = 8)
hist(unlist(normCounts), breaks=100000, las=1,
     xlab="normalized counts", ylab="Occurrences", col="grey",
     main=paste(parameters$recountID, "Normalised count distrib"),
     xlim=c(0,q0.95))
abline(v=mean(unlist(normCounts)), lwd=1, col="darkgreen") # mark the mean
abline(v=median(unlist(normCounts)), lwd=1, col="blue") # mark the median
abline(v=mean(unlist(trimmed)), lwd=1, col="purple") # mark the trimmed mean
legend("topright",lwd=2,
       legend=c("mean", "median", "trimmed mean"),
       col=c("darkgreen", "blue", "purple"))
silence <- dev.off() ## Note: I send the result to a variable to avoid verbosity

## Interpretation:  a huge difference between the mean and median count values.
## This means that the distribution is highly right-skewed.
## In particular, the mean is pulled on the right by the outliers.
## This effect is slighly reduced when we discard the zeros and the values
## above percentile 965, but there is still a skew.
##
## The log2 transformation reduces this problem. Thenw e always work with log2 counts.


## Histogram of log2_transformed normalized counts
pdf(file=file.path(
  dir.log2Impact,
  paste(parameters$recountID, "_log2counts_norm-perc75_hist", ".pdf", sep = "")),
  width = 7, height = 5)
hist(unlist(log2norm$Counts), breaks=100, las=1,
     xlab="log2(counts)", ylab="Occurrences", col="grey",
     main=paste(parameters$recountID, "Normalised log2-count distrib"))
abline(v=mean(unlist(log2norm$Counts)), lwd=2, col="darkgreen") # mark the median
abline(v=median(unlist(log2norm$Counts)), lwd=2, col="blue") # mark the median
legend("topright",lwd=2,
       legend=c("mean", "median"),
       col=c("darkgreen", "blue"))
#dev.copy(png , "Normalised_log2-count_distribution")
silence <- dev.off()

## JvH: Mustafa: how do you interpret the difference between
## the mean and the median for the raw counts ?
pdf(file=file.path(
  dir.log2Impact,
  paste(parameters$recountID, "_log2counts_norm-perc75_mean_vs_median", ".pdf", sep = "")),
  width = 8, height = 8)
plot(apply(log2norm$Counts, 1, mean),
     apply(log2norm$Counts, 1, median),
     main="log2norm$Counts: Median versus mean",
     xlab="Mean log2(norm counts) per sample",
     ylab="Median log2(norm counts) per sample",
     col=sampleColors,
     las=1,
     panel.first=grid())
legend("bottomright", legend=names(classColors), col=classColors, pch=1, cex=0.8)
silence <- dev.off()

## JvH: Mustafa: how do you interpret the difference between
## the mean and the third quartile ?
pdf(file=file.path(
  dir.figures,
  paste(parameters$recountID, "_counts_mean_vs_Q3", ".pdf", sep = "")),
  width = 8, height = 8)
plot(x=apply(rawCounts1, 1, mean),
     y=signif(digits=3, apply(rawCounts1, 1, quantile, 0.75)),
     main="raw counts: Percentile 75 versus mean",
     xlab="Mean counts per sample",
     ylab="Percentile 75",
     col=sampleColors,
     las=1,
     panel.first=grid())
legend("topleft", legend=names(classColors), col=classColors, pch=1, cex=0.8)
silence <- dev.off()

# ## Plot a summary of the most popular normalising factor
# pdf(file=file.path(
#   dir.figures,
#   paste(parameters$recountID, "_log2counts_norm-perc75_scaling_factors", ".pdf", sep = "")),
#   width = 12, height = 12)
# plot(log2norm$sampleStats[, c("mean", "median", "Q3")],
#      main="Choice of the scaling factor",
#      col=sampleColors,
#      las=1)
# silence <- dev.off()

## JvH: Mustafa, please pay attention to this table of sample-wise statistics.
## Analyse the first quantile and the median. What do you conclude ? Do you think
## that these could be suitably used as normalisation factors ?
# View(log2norm$sampleStats)

## Draw a box plot of the log2-transformed normalised counts.
## With normalisation method quantile=0.75, the third quantiles should be
## perfectly aligned.
## NOTE: THIS MAKES R GET STUCK -> I COMMENT IT.
# par(mfrow=c(1,2))
# boxplot(t(normCounts), horizontal=TRUE,  cex.axis=0.4, las=1,
#         col=sampleColors,
#         main="p75-normalised counts per sample")
# boxplot(t(log2norm$Counts), horizontal=TRUE,  cex.axis=0.4, las=1,
#         col=sampleColors,
#         main="log2-transformed p75--normalised counts per sample")
# par(mfrow=c(1,1))


################################################################
## Miscalssification error rates without any traetment onto the
## raw count.

################################################################
## Iterate to measure the testing error rate for the KNN classifier.

## Iterate over classifiers
if (parameters$compute) {
  train.test.results <- list()

  message.with.time("\n\nTrain/test computations with all variables")

  for (classifier in parameters$classifiers) {
    if (classifier == "lda") {
      message.with.time("Skipping LDA for all variables")
    } else {
      # message.with.time("Classifier: ", classifier)

      ## Iterate over data types
      for (data.type in parameters$data.types) {
        # message.with.time("Data type: ", data.type)

        ## choose counts according to data type
        if (data.type == "raw") {
          counts <- rawCounts
        } else if (data.type == "norm") {
          counts <- normCounts
        } else if (data.type == "log2") {
          counts <- log2norm$Counts
        } else {
          stop(data.type, " is not a valid data type. Supported: raw, norm, log2")
        }

        ## Iterate over permutation status
        for (permute in c(FALSE, TRUE)) {
          # message.with.time("Permute class labels: ", permute)

          ## Define experiment prefix
          exp.prefix <- paste(sep="_", classifier, parameters$recountID, data.type, "all")
          if (permute) {
            exp.prefix <- paste(sep="_", exp.prefix, perm.prefix)
          }
          message (format(Sys.time(), "%Y-%m-%d_%H%M%S"), "\t", "Experiment prefix: ", exp.prefix)

          train.test.results[[exp.prefix]] <-
            one.experiment(countTable=counts,
                           data.type=data.type,
                           classifier=classifier,
                           classes=classes,
                           variable.type = "all",
                           permute = permute,
                           trainingProportion = parameters$trainingProportion,
                           file.prefix = exp.prefix)
        } ## End iterate over permutation status
      } ## End iterate over data types
    } ## End iterate over classifiers
  }
} else {
  message.with.time("Skipping train/test computations with all variables")
}


################################################################
####   Variable ordering based on differential expression   ####
################################################################
if (parameters$compute) {
  ## Run differential analysis with DESeq2 and edgeR to define variable order
  message.with.time("Running DESeq2 and edgeR to define variable ordering")
  DEG.DESeq2 <- DEGordering(loaded$countTable, loaded$classes, method = "DESeq2")
  DEG.edgeR  <- DEGordering(loaded$countTable, loaded$classes, method = "edgeR")

  sorted.log2.transformed.edgeR <- loaded$countTable[, DEG.DESeq2$geneOrder]
  sorted.log2.transformed.DESeq2 <- loaded$countTable[, DEG.edgeR$geneOrder]
} else {
  message.with.time("Skipping DEG detection with edgeR and DESeq2. ")
}

################################################################
## Second experiment: measure hit rates with increasing number of variables ordered by DEG p-value
################################################################
if (parameters$compute) {
  # iteratedKnn.test.edgeR.allVariables <- data.frame()
  # iteratedKnn.test.DESeq2.allVariables <- data.frame()

  for (classifier in parameters$classifiers) {
    for (permute in c(FALSE, TRUE)) {
      for (deg.method in parameters$deg.methods) {
        DEG <- get(paste(sep="", "DEG.", deg.method))

        v  <- 1
        for(v in 1:length(parameters$nb.variables)){
          varnb <- parameters$nb.variables[v]

          ## For the time being we do this experiment only with log2 normalised counts
          ## since we saw that it improves the result with all variables
          data.type <- "log2"
          counts <- (log2norm$Counts[,DEG$geneOrder[1:varnb]])
          ## dim(counts)

          ## Define experiment prefix
          variable.type <- paste(sep="_", "DEG", deg.method, "top", varnb)
          exp.prefix <- paste(sep="_", classifier, parameters$recountID, data.type, variable.type)
          if (permute) {
            exp.prefix <- paste(sep="_", exp.prefix, perm.prefix)
          }
          message (format(Sys.time(), "%Y-%m-%d_%H%M%S"), "\t", "Experiment prefix: ", exp.prefix)

          train.test.results[[exp.prefix]] <-
            one.experiment(countTable=counts,
                           data.type=data.type,
                           classifier=classifier,
                           classes=classes,
                           variable.type = variable.type,
                           permute = permute,
                           file.prefix = exp.prefix)
        }
      }
    }
  }
}  else {
  message.with.time("Skipping train/test analysis with DEG-ordered top variables")
}

## For each experiment, the results were stored in an element of the list train.test.results.
## Test the list of experiments done.
names(train.test.results)

################################################################
## Gather error rates for experiments with all variables and draw
## a box plot.
all.variables.testing.error.rates <- data.frame()
all.variables.experiments <- grep (pattern = "_all", x = names(train.test.results), value = TRUE)
if (is.null(all.variables.experiments) || (length(all.variables.experiments)<1)) {
  message.with.time("Not a single experiment for all variables composition")
} else {


  for (exp in 1:length(all.variables.experiments)) {
    experiment <- all.variables.experiments[exp]
    exp.result <- train.test.results[[experiment]]
    #  all.variables.testing.error.rates[,experiment] <- exp.result$testing.error.rate
    if (exp == 1) {
      all.variables.testing.error.rates <- exp.result$testing.error.rate
    } else {
      all.variables.testing.error.rates <- cbind(all.variables.testing.error.rates, exp.result$testing.error.rate)
    }
  }
  colnames(all.variables.testing.error.rates) <- all.variables.experiments
  rownames(all.variables.testing.error.rates) <- 1:nrow(all.variables.testing.error.rates)

  ## Define file to store the boxplot
  file.prefix <- "all_variables"
  boxplot.file <- file.path(dir.figures, paste(sep="", file.prefix, "_boxplot.pdf"))
  message.with.time("Boxplot file: ", boxplot.file)
  pdf(file=boxplot.file, width = 8, height = 10)

  ## Define parameters for the boxplot
  # pdf(file=boxplot.file, width = 7, height = 5)
  save.margins <- par("mar")
  par(mar = c(10, 5, 7, 1))
  testTable.colors <- c(1:ncol(all.variables.testing.error.rates))
  names(testTable.colors) <- colnames(all.variables.testing.error.rates)
  testTable.colors[grep(names(testTable.colors),pattern = "_raw")] <- "pink" ## Gray for permutation tests
  testTable.colors[grep(names(testTable.colors),pattern = "_norm")] <- "yellow" ## Gray for permutation tests
  testTable.colors[grep(names(testTable.colors),pattern = "_log2")] <- "green" ## Gray for permutation tests
  testTable.colors[grep(names(testTable.colors),pattern = "_permLabels")] <- "gray" ## Gray for permutation tests
  boxplot(all.variables.testing.error.rates,
          horizontal = FALSE, las= 2, cex.axis = 0.8, col= testTable.colors,
          ylab= "Misclassification rate",
          ylim=c(0,1),
          main= paste("Misclassification rates", "\nAll variables,", parameters$iterations ,"iterations"  ))
  abline(h=seq(from=0,to=1,by=0.1), col="darkgrey", lty="dotted")
  par(mar=save.margins)
  silence <- dev.off()
}

################################################################
##
## Gather error rates for experiments with ordered variables
## and generate boxplots to show the effect of variable number.
## Do it separately for each DEG method and each classifier.
##
################################################################
for (classifier in parameters$classifiers) {
  for (deg.method in parameters$deg.methods) {
    for (permute in c(FALSE, TRUE)) {
      data.type <- "log2"
      variable.type <- "DEG-ordered"

      main <- paste(sep="", classifier, "; ", deg.method, " ", variable.type, " variables; ", data.type, " counts; ","\n", parameters$iterations, " iterations")
      if (permute) {
        main <- paste(main, "; Permuted class labels")
      }
      message.with.time("Generating error rate box plot: ", main)

      ## Gather error rates for experiments with all variables
      selected.experiments <- paste(sep="_", classifier, data.type, "DEG", deg.method, "top", parameters$nb.variables)
      if (classifier == "lda") {
        message.with.time("Skipping all variables for LDA")
      } else {
        selected.experiments <- c(selected.experiments, paste(sep="_", classifier, data.type, "all"))
      }
      if (permute) {
        selected.experiments <- paste(sep="_", selected.experiments, perm.prefix)
      }

      ## Collect the error rates for the selected experiments
      top.variables.testing.error.rates <- data.frame()
      for (exp in 1:length(selected.experiments)) {
        experiment <- selected.experiments[exp]
        exp.result <- train.test.results[[experiment]]
        if (exp == 1) {
          top.variables.testing.error.rates <- data.frame(exp.result$testing.error.rate)
          colnames(top.variables.testing.error.rates) <- experiment
        } else {
          top.variables.testing.error.rates[, experiment] <- exp.result$testing.error.rate
        }
      }

      ## Define file to store the boxplot
      file.prefix <- paste(sep = "_", classifier, deg.method, "ordered")
      if (permute) {
        file.prefix <- paste(sep="_", file.prefix, perm.prefix)
      }
      boxplot.file <- file.path(dir.figures, paste(sep="", file.prefix, "_DEG-ordered_top-variables_boxplot.pdf"))
      message.with.time("Boxplot file: ", boxplot.file)
      pdf(file=boxplot.file, width = 10, height = 8)

      ## Define parameters for the boxplot
      save.margins <- par("mar")
      par(mar = c(6, 5, 6, 1))
      boxplot.names <- parameters$nb.variables
      boxplot.names <- paste("top", boxplot.names)
      if (classifier != "lda") {
        boxplot.names <- append(boxplot.names, paste(sep="", "all (", ncol(rawCounts), ")"))
      }


      testTable.colors <- c(1:ncol(top.variables.testing.error.rates))
      names(testTable.colors) <- colnames(top.variables.testing.error.rates)
      testTable.colors[grep(names(testTable.colors),pattern = "_raw")] <- "pink" ## Gray for permutation tests
      testTable.colors[grep(names(testTable.colors),pattern = "_norm")] <- "yellow" ## Gray for permutation tests
      testTable.colors[grep(names(testTable.colors),pattern = "_log2")] <- "green" ## Gray for permutation tests
      testTable.colors[grep(names(testTable.colors),pattern = "_permLabels")] <- "gray" ## Gray for permutation tests
      boxplot(top.variables.testing.error.rates,
              horizontal = FALSE, las= 2, cex.axis = 1, col= testTable.colors,
              ylab= "Misclassification rate",
              names = boxplot.names,
              ylim=c(0,1),
              main= main)
      abline(h=seq(from=0,to=1,by=0.1), col="darkgrey", lty="dotted")
      par(mar=save.margins)
      silence <- dev.off()
    }
  }
}



######################################################################
## Plot contingency tables
#####################################################################

## TO BE DONE

######################################################################
## Clustering of genes and samples (based on significant DEG)
#####################################################################

## TO BE DONE


######################################################################
## Functional enrichment analysis
#####################################################################

## TO BE DONE


stop("Mustafa, the next analyses should now be done with the same one.experiment() function")


###########################################################################################################
########## Testing KNN  & RF classifier performance with the DESeq2 and edgeR ordaring variables for raw count
###########################################################################################################
## we thus will using the ordered variables which are the most significance resulting from DEG,
## concering with padj.

second.experiment <- function( countTable = rawCounts,
                               classes = classes,
                               classifier = "knn",
                               data.type = "raw",
                               variable.type = "top_20000",
                               permute = F,
                               deg.method = "DESeq2"
)
{
  message.with.time(classifier, " :classifier; ", deg.method , " :DEG method; ", data.type,":counts; ",
                    variable.type, ":variable number" )

  if(deg.method == "DESeq2") {
    DEG.DESeq <- DEGordering(countTable , classes, method = "DESeq2")
    sorted.countTable.DESeq <- countTable[, DEG.DESeq$geneOrder]
  }else {
    DEG.edgeR  <- DEGordering(countTable , classes, method = "edgeR")
    sorted.countTable.edgeR <- countTable[,  DEG.edgeR$geneOrder]
  } # end of the IF for checking the DEG( Deffirential expression genes)


  v  <- 1
  for(v in 1:length(parameters$nb.variables)){

    varnb <- parameters$nb.variables[v]
    if( deg.method == "DESeq2"){
      message.with.time("DESeq2-ordered ", v, "/", length(nb.variables), ": ",varnb," variables")
      testTable <- data.frame()
      testTable.by.varnb <- list()
      var.label <- paste(sep="_", "DESeq2_top", varnb)

      for(i in 1:parameters$iterations){
        test.DESeq <-MisclassificationEstimate(countTable = sorted.countTable.DESeq[, 1:varnb], classes = classes, classifier = classifier)
        testTable <- rbind(testTable, test.DESeq$stats)
      } # end for iteration

      # Store the full result table for each var number
      testTable.by.varnb[[var.label]] <- testTable

      # Collect the error rates in a single table with one row per iteration and one col per number of variables
      if (v == 1) {

        testTable.error.rate <- data.frame(testTable$testing.error.rate)
        colnames(testTable.error.rate) <- var.label
      } else {
        testTable.error.rate[, var.label] <-
          testTable$testing.error.rate
      } # end of else V==1


    }else { ## if the DEG is edgeR
      message.with.time("edgeR-ordered ", v, "/", length(parameters$nb.variables), ": ",varnb," variables")
      testTable.by.varnb <- list()
      var.label <- paste(sep="_", "edgeR_top", varnb)

      for(i in 1:parameters$iterations){
        test.edgeR <-MisclassificationEstimate(countTable = sorted.countTable.edgeR[, 1:varnb], classes = classes ,classifier = classifier )
        testTable <- rbind(testTable, test.edgeR$stats)
      } # end for iteration

      # Store the full result table for each var number
      testTable.by.varnb[[var.label]] <- testTable

      # Collect the error rates in a single table with one row per iteration and one col per number of variables
      if(v == 1){
        testTable.error.rate <- data.frame( testTable$testing.error.rate)
        colnames(testTable.error.rate) <- var.label
      }else{
        testTable.error.rate[,var.label] <- testTable$testing.error.rate
      }# end of else for V== 1

    }# end if DEG is edgeR

  }# end of for variables vector

  testTable.error.rate[,paste("All variables (", ncol(countTable), ")")] <-
    data.frame(testTable$testing.error.rate)

  testTable.error.rate[, paste("All variables (",ncol(countTable), ")")] <-
    data.frame( testTable$testing.error.rate)

  ## extracting the results for the testing and learning error rate in tables and figures
  testTable.prefix <- paste(sep = "_", classifier, data.type, variable.type, deg.method)
  if(classifier == "knn"){
    testTable.file <- file.path(table.dirs[classifier], paste(sep="", testTable.prefix, ".tsv"))
    errorSummary.file <- file.path(table.dirs[classifier], paste(sep="", testTable.prefix, "_summary.tsv"))
    boxplot.file <- file.path(
      figure.dirs[classifier],
      paste(sep="", testTable.prefix, "_R", parameters$iterations, "_learning_vs_test_error_boxplot.pdf"))
  }else {
    testTable.file <- file.path(table.dirs[classifier], paste(sep="", testTable.prefix, ".tsv"))
    errorSummary.file <- file.path(table.dirs[classifier], paste(sep="", testTable.prefix, "_summary.tsv"))
    boxplot.file <- file.path(
      figure.dirs[classifier],
      paste(sep = "", testTable.prefix , "_R", parameters$iterations, "_learning_vs_test_error_boxplot.pdf"))
  }

  if (parameters$save.tables) {
    write.table(file=testTable.file,
                x = testTable.by.varnb,
                row.names = TRUE, col.names = NA, sep="\t", quote=FALSE)
  } else {
    message.with.time("Skipping saving of test table")
  }

  if (parameters$save.tables) {
    ## Export a summary about the error rate
    write.table(tidy( summary(testTable$testing.error.rate)),
                quote=FALSE, sep="\t", row.names=TRUE, col.names = NA,
                file = errorSummary.file)
  } else {
    message.with.time("Skipping saving of test summary table")
  }

  ## Plot a box plot of training versus testing error
  pdf(file=boxplot.file, width = 7, height = 5)
  save.margins <- par("mar")
  par(mar = c(3, 10, 7, 1))
  testTable.colors <- c(1:length(testTable.error.rate))
  names(testTable.colors)<- names(testTable.error.rate)
  testTable.colors[grep(names(testTable.colors),pattern = "top")] <- "#9900FF"
  testTable.colors[grep(names(testTable.colors),pattern = "All")] <- "gray"
  main <- paste(sep="", classifier, "; ", deg.method, "; " , data.type, " counts; ", variable.type, " variables")
  if (permute) {
    main <- paste(sep="", main, "; permuted labels")
  }

  if (variable.type == "top") {
    boxplot(testTable.error.rate , horizontal = T ,las= 1 ,cex.axis = 0.8,col= testTable.colors,
            xlab= "Misclassification Rate ",
            ylim=c(0,1),
            main= paste("Misclassification Rate\n", main , "\n", parameters$iterations ,"iteration"  ))
    silence <- dev.off()
  } else {
    boxplot(testTable[, c("training.error.rate", "testing.error.rate")],
            ylim=c(0,1),
            main=main)
    silence <- dev.off()
  } # end of else for boxplot
  return(testTable)

}# end function second.experiment

## Run KNN classifier with all variables
knnDEGrawTestTable <- second.experiment(countTable=rawCounts,
                                        data.type="raw",
                                        classifier="knn",
                                        deg.method = "DESeq2",
                                        classes=classes,
                                        variable.type = "all",
                                        permute = FALSE)
## Run KNN classifier with 20000 top variables
knnDEGrawTestTable <- second.experiment(countTable=rawCounts,
                                        data.type="raw",
                                        classifier="knn",
                                        deg.method = "DESeq2",
                                        classes=classes,
                                        variable.type = "top20000",
                                        permute = FALSE)

## Run RF classifier with all variables
knnDEGrawTestTable <- second.experiment(countTable=rawCounts,
                                        data.type="raw",
                                        classifier="rf",
                                        deg.method = "DESeq2",
                                        classes=classes,
                                        variable.type = "all",
                                        permute = FALSE)
## Run RF classifier with top variables
knnDEGrawTestTable <- second.experiment(countTable=rawCounts,
                                        data.type="raw",
                                        classifier="rf",
                                        deg.method = "DESeq2",
                                        classes=classes,
                                        variable.type = "top20000",
                                        permute = FALSE)


#########################################################################################################
if (parameters$compute) {

  ## Instantiate result tables and lists
  KNNlog2TestTable.edgeRordering <- list()
  edgeR.error.table <- data.frame()

  KNNlog2TestTable.DESeq2ordering <- list()
  DESeq2.error.table <- data.frame()

  sorted.log2.transformed.edgeR <- log2norm$Counts[, DEG.DESeq2$geneOrder]
  sorted.log2.transformed.DESeq2 <- log2norm$Counts[, DEG.edgeR$geneOrder]

  iteratedKnn.test.edgeR.allVariables <- data.frame()
  iteratedKnn.test.DESeq2.allVariables <- data.frame()

  v <- 1 ## for quick test
  for( v in 1:length(parameters$nb.variables)){
    varnb <- parameters$nb.variables[v]
    message.with.time("classifier with ordered Variables ", v, "/", length(parameters$nb.variables) , " variables")

    # for(i in 1:parameters$iterations){
    #   message.with.time("RF classifier with ", i, "/", parameters$iterations , " iterations")
    #   test.edgeR <-MisclassificationEstimate(countTable = sorted.log2.transformed.edgeR[, 1:varnb], classes = classes ,classifier = "knn" )
    #   test.DESeq2 <-MisclassificationEstimate(countTable = sorted.log2.transformed.DESeq2[, 1:varnb], classes = classes, classifier = "knn")
    #   iteratedKnn.test.edgeR.allVariables <- rbind(iteratedKnn.test.edgeR.allVariables, test.edgeR$stats)
    #   iteratedKnn.test.DESeq2.allVariables <- rbind(iteratedKnn.test.DESeq2.allVariables, test.DESeq2$stats)
    #
    # } # end for iteration

    # Run KNN classifier with all variables, log2-transformed normalised with percentile 75
    KNNlog2TestTable.edgeRordering[[varnb]] <-  one.experiment(countTable=sorted.log2.transformed.edgeR[, 1:varnb],
                                                               data.type="log2normP75",
                                                               classifier=classifier,
                                                               classes=classes,
                                                               variable.type = paste(sep = '', 'top', varnb),
                                                               permute = FALSE)

    # meanOfEachIteration  <- apply( iteratedKnn.test.edgeR.allVariables["testing.error.rate"], 2, mean)
    # iteratedmeanOfEachIterationAll <- append( iteratedmeanOfEachIteration ,meanOfEachIteration )
  }
} else {
  message.with.time("Skipping testing KNN classifier with ordered variables for real dataset")
}



###########################################################################################################
########## Testing KNN classifier performance with the DESeq2 and edgeR ordaring variables for real data set
###########################################################################################################
## we thus will using the ordered variables which are the most significance resulting from DEG,
## concering with padj.

if(parameters$compute){
  message.with.time("KNN classifier with DESeq2 and edgeR ordaring real data set , ", parameters$iterations, " iterations.")

  DEG.DESeq2 <- DEGordering(loaded$countTable , loaded$classes, method = "DESeq2")
  DEG.edgeR  <- DEGordering(loaded$countTable, loaded$classes, method = "edgeR")

  sorted.log2.transformed.edgeR <- loaded$countTable[, DEG.DESeq2$geneOrder]
  sorted.log2.transformed.DESeq2 <- loaded$countTable[, DEG.edgeR$geneOrder]

  iteratedKnn.test.edgeR.allVariables <- data.frame()
  iteratedKnn.test.DESeq2.allVariables <- data.frame()

  edgeR.error.table <- data.frame()
  DESeq2.error.table <- data.frame()


  for( v in 1:length(parameters$nb.variables)){
    message.with.time("RF classifier with ordered Variables", v, "/", length(parameters$nb.variables) , "variables")
    varnb <- parameters$nb.variables[v]

    for(i in 1:parameters$iterations){
      message.with.time("RF classifier with", i, "/", parameters$iterations , "iteration")
      test.edgeR <-MisclassificationEstimate(countTable = sorted.log2.transformed.edgeR[, 1:varnb], classes = classes ,classifier = "knn" )
      test.DESeq2 <-MisclassificationEstimate(countTable = sorted.log2.transformed.DESeq2[, 1:varnb], classes = classes, classifier = "knn")
      iteratedKnn.test.edgeR.allVariables <- rbind(iteratedKnn.test.edgeR.allVariables, test.edgeR$stats)
      iteratedKnn.test.DESeq2.allVariables <- rbind(iteratedKnn.test.DESeq2.allVariables, test.DESeq2$stats)


    } # end for iteration
    meanOfEachIteration  <- apply( iteratedKnn.test.edgeR.allVariables["testing.error.rate"], 2, mean)
    iteratedmeanOfEachIterationAll <- append( iteratedmeanOfEachIteration ,meanOfEachIteration )
  } # end for number of variables
}else{
  message.with.time("Skipping testing KNN classifier with ordered variables for real dataset")
}

## Hereafter we noticed that ordaring of variable has achieved much improvement in the calssification
## performing the Differential expression analysis between our two classes (Heparinised.blood & Bone marrow)
## will foster the efficiancy's classifier

## results after we applying edgeR DE and arraning the variables with the most significance to study
## the effect of such ordering onto the efficaincy of the KNN classifier.
#iterated.test.edgeR

### Results of KNN real data set with edgeR ordaring variables
if (parameters$save.tables) {

  write.table(file = file.path(table.dirs[classifier], "testing_KNN_with_edgeR_ordered_variables_for_real_dataset.tsv"),
              quote=FALSE, sep="\t", row.names=TRUE, col.names = NA,
              x = iteratedKnn.test.edgeR)
}

summaryKnnedgeRTidy <- tidy(summary(iteratedKnn.test.edgeR$testing.error.rate))
if (parameters$save.tables) {
  write.table(file = file.path(table.dirs[classifier], "KNN_Summary_with_edgeR_ordered_valriable_for_real_dataset.tsv"),
              quote=FALSE, sep="\t", row.names=TRUE, col.names = NA,
              x = summaryKnnedgeRTidy)
}


pdf(file = file.path(figure.dirs[classifier], paste("KNN_testing_with_edgeR_ordered_variables_for_real_dataset", ".pdf", sep = "")),
    width =9 , height = 5 )
boxplot(iteratedKnn.test.edgeR[,c("training.error.rate", "testing.error.rate")],
        ylim = c(0,1),
        main= "Testing KNN classifier with edgeR ordered variables for real dataset")
silence <- dev.off()

########################################################
## results after we applying DESeq2 DE and arraning the variables with the most significance to study
## the effect of such ordering onto the efficaincy of the KNN classifier.

#### Results of KNN real data set with DESseq2 ordaring variables
if (parameters$save.tables) {
  write.table(file = file.path(table.dirs[classifier], "Testing_KNN_classifier_DESeq2_ordered_variables_for_real_dataset.tsv"),
              quote=FALSE, sep="\t", row.names=TRUE, col.names = NA,
              x = iteratedKnn.test.DESeq2)
}

summaryKnnDESeq2Tidy <- tidy(summary(iteratedKnn.test.DESeq2$testing.error.rate))
if (parameters$save.tables) {
  write.table(file = file.path(figure.dirs[classifier], "testing_KNN_with_DESeq2_ordered_variables_for_real_dataset.tsv"),
              quote=FALSE, sep="\t", row.names=TRUE, col.names = NA,
              x = summaryKnnDESeq2Tidy)
}

pdf(file = file.path(figure.dirs[classifier], paste( "KNN_testing_with_DESeq2_ordered_for_real_dataset", ".pdf", sep = "")),
    width =9 , height = 5 )
boxplot(iteratedKnn.test.DESeq2[,c("training.error.rate", "testing.error.rate")],
        ylim = c(0,1),
        main= "Testing KNN classifier with DESeq2 ordered variable for real dataset")
silence <- dev.off()


###########################################################################################################
##########  Testing KNN classifier performance with the DESeq2 and edgeR ordaring variables for
##########  Log2-transformed-Normalized count dataset
###########################################################################################################
## we thus will using the ordered variables which are the most significance resulting from DEG,
## concering with padj.
if(parameters$compute) {
  message.with.time("KNN classifier with DESeq2 and edgeR ordaring real data set , ", parameters$iterations, " iterations.")

  DEG.DESeq2 <-  DEGordering(loaded$countTable , loaded$classes, method = "DESeq2")
  DEG.edgeR  <- DEGordering(loaded$countTable, loaded$classes, method = "edgeR")

  sorted.log2.transformed.edgeR <-log2norm$Counts[, DEG.DESeq2$geneOrder]
  sorted.log2.transformed.DESeq2 <- log2norm$Counts[, DEG.edgeR$geneOrder]

  iteratedKnn.test.edgeR.Log2Norm <- data.frame()
  iteratedKnn.test.DESeq2.Log2Norm <- data.frame()

  edgeR.error.table <- data.frame()
  DESeq2.error.table <- data.frame()
  iteratedVarnb <- vector()
  iteratedmeanOfEachIteration <- vector()
  for(i in 1:length(parameters$nb.variables)){
    varnb <- parameters$nb.variables[i]

    for(i in 1:parameters$iterations){
      test.edgeR <-MisclassificationEstimate(countTable = sorted.log2.transformed.edgeR[, 1:varnb], classes = classes ,classifier = "knn" )
      test.DESeq2 <-MisclassificationEstimate(countTable = sorted.log2.transformed.DESeq2[, 1:varnb], classes = classes, classifier = "knn")
      iteratedKnn.test.edgeR.Log2Norm <- rbind(iteratedKnn.test.edgeR.Log2Norm, test.edgeR$stats)
      iteratedKnn.test.DESeq2.Log2Norm <- rbind(iteratedKnn.test.DESeq2.Log2Norm, test.DESeq2$stats)


    } # end for iteration
    meanOfEachIteration  <- apply( iteratedKnn.test.edgeR.Log2Norm["testing.error.rate"], 2, mean)
    iteratedmeanOfEachIteration <- append( iteratedmeanOfEachIteration ,meanOfEachIteration )
    iteratedVarnb <- append(iteratedVarnb, varnb )
  } # end for number of variables
} else {
  message.with.time("Skipping testing KNN classifier with log2-normalised with DESeq2 and edgeR ordering variables for real dataset ")
}
## Hereafter we noticed that ordaring of variable has achieved much improvement in the calssification
## performing the Differential expression analysis between our two classes (Heparinised.blood & Bone marrow)
## will foster the efficiancy's classifier

## results after we applying edgeR DE and arraning the variables with the most significance to study
## the effect of such ordering onto the efficaincy of the KNN classifier.
#iterated.test.edgeR



##################################################################################################
########## Testing KNN classifier performance with the DESeq2 and edgeR ordaring permuted variables of dataset
##################################################################################################

################################################################
## Do a negative control: estimate the misclassification rate
## When we foold the classifier by giving it randomly permuted classes.
## This gives us an estimation of the misclassification rate expected by chance.
## The question is "how much better does the classifier perform with real classes
## than with randomized classes" ?
if(parameters$compute) {
  message.with.time("analysis KNN classifier with DESeq2 and edgeR permuted ordaring data set , ", parameters$iterations, " iterations.")

  iteratedKnn.test.edgeRPermuted <- data.frame()
  iteratedKnn.test.DESeq2Permuted <- data.frame()

  edgeR.error.table <- data.frame()
  DESeq2.error.table <- data.frame()
  for(i in 1:length(parameters$nb.variables)){
    varnb  <-parameters$nb.variables[i]

    for(i in 1:parameters$iterations){

      test.edgeR <- MisclassificationEstimate(countTable = sorted.log2.transformed.edgeR[, 1:varnb], classes = sample(classes) ,classifier = "knn" )
      test.DESeq2 <- MisclassificationEstimate(countTable = sorted.log2.transformed.DESeq2[, 1:varnb], classes = sample(classes), classifier = "knn")

      iteratedKnn.test.edgeRPermuted  <- rbind(iteratedKnn.test.edgeRPermuted  , test.edgeR$stats)
      iteratedKnn.test.DESeq2Permuted <- rbind(iteratedKnn.test.DESeq2Permuted , test.DESeq2$stats)

    } # end for iteration
  }   # end for number of variables
} else {
  message.with.time("Testing KNN classifier with DESeq2 & edgeR ordered variables for permuted dataset")
}
## these are studing the effect of permuted classes where we tried fooled classifier by give it
## randomly permuted classe to know if this will effert the efficacy for classifier in calssification
## process

### Results of KNN permuted with edgeR ordaring variables
write.table(file = file.path(table.dirs[classifier], "KNN_testing_with_edgeR_ordered_variables_for_permuted_dataset.tsv"),
            quote=FALSE, sep="\t", row.names=TRUE, col.names = NA,
            x = iteratedKnn.test.edgeRPermuted )
summaryedgeRPermutedTidy <- tidy(summary(iteratedKnn.test.edgeRPermuted$testing.error.rate))
write.table(file = file.path(figure.dirs[classifier], "KNN_testing_edgeR_with_ordered_variables_for_permuted_dataset.tsv"),
            quote=FALSE, sep="\t", row.names=TRUE, col.names = NA,
            x = summaryedgeRPermutedTidy)
pdf(file = file.path(table.dirs[classifier], paste( "KNN_testing_with_edgeR_ordered_variables_for_permuted_dataset", parameters$iterations,
                                           ".pdf")), width = 9 , height = 5)
boxplot(iteratedKnn.test.edgeRPermuted[, c("training.error.rate", "testing.error.rate")],
        ylim= c(0,1), main="Testing KNN classifier with edgeR ordered variables for permuted dataset")
silence <- dev.off()

### Results of KNN permuted with DESeq2 ordaring variables
write.table(file = file.path(table.dirs[classifier], "KNN_testing_with_DESeq2_ordered_variables_for_permuted_dataset.tsv"),
            quote=FALSE, sep="\t", row.names=TRUE, col.names = NA,
            x = iteratedKnn.test.DESeq2Permuted)
summaryDESeq2PermutedTidy <- tidy(summary(iteratedKnn.test.DESeq2Permuted$testing.error.rate))
write.table(file = file.path(figure.dirs[classifier], "KNN_testing_DESeq2_ordered_varailes_for_permuted_dataset.tsv"),
            quote=FALSE, sep="\t", row.names=TRUE, col.names = NA,
            x = summaryDESeq2PermutedTidy)
pdf(file = file.path(dir.figures, paste( "KNN_testing_DESeq2_ordered_varibles_for_permuted_real_dataset", parameters$iterations,
                                         ".pdf")), width = 9 , height = 5)
boxplot(iteratedKnn.test.DESeq2Permuted[, c("training.error.rate", "testing.error.rate")],
        ylim= c(0,1), main="Testing KNN Classifier with DESeq2 ordered variables for permuted dataset")
silence <- dev.off()


###########################################################################################################
###### Testing random Forest classifier performance with the DESeq2 and edgeR ordaring variables for #####
###### real data set                                                                                 ######
###########################################################################################################
## we thus will using the ordered variables which are the most significance resulting from DEG,
## concering with padj.
if(parameters$compute){
  message.with.time("analysis RF classifier with DESeq2 and edgeR ordaring real data set , ", parameters$iterations, " iterations.")

  # DEG.DESeq2 <- DEGordering(loaded$countTable , loaded$classes, method = "DESeq2")
  # DEG.edgeR  <- DEGordering(loaded$countTable, loaded$classes, method = "edgeR")

  sorted.log2.transformed.edgeR <-log2norm$Counts[, DEG.DESeq2$geneOrder]
  sorted.log2.transformed.DESeq2 <- log2norm$Counts[, DEG.edgeR$geneOrder]

  iteratedRf.test.edgeR <- data.frame()
  iteratedRf.test.DESeq2 <- data.frame()

  edgeR.error.table <- data.frame()
  DESeq2.error.table <- data.frame()


  for(i in 1:length(parameters$nb.variables)){
    varnb <- parameters$nb.variables[i]

    for(i in 1:parameters$iterations){
      test.edgeR <-MisclassificationEstimate(countTable = sorted.log2.transformed.edgeR[, 1:varnb], classes = classes ,classifier = "rf" )
      test.DESeq2 <-MisclassificationEstimate(countTable = sorted.log2.transformed.DESeq2[, 1:varnb], classes = classes, classifier = "rf")
      iteratedRf.test.edgeR <- rbind(iteratedRf.test.edgeR, test.edgeR$stats)
      iteratedRf.test.DESeq2 <- rbind(iteratedRf.test.DESeq2, test.DESeq2$stats)

    } # end for iteration
  } # end for number of variables
}else {
  message.with.time("Skipping RF classifier testing with log2-normalised for real dataset")
}
## Hereafter we noticed that ordaring of variable has achieved much improvement in the calssification
## performing the Differential expression analysis between our two classes (Heparinised.blood & Bone marrow)
## will foster the efficiancy's classifier

## results after we applying edgeR DE and arranging the variables with the most significance to study
## the effect of such ordering onto the efficaincy of the RF classifier.
#iterated.test.edgeR

### Results of FR real data set with edgeR ordaring variables
write.table(file = file.path(table.dirs[classifier], "Testing_RF_with_log2-nrmalised_edgeR_ordered_variables_for_real_dataset.tsv"),
            quote=FALSE, sep="\t", row.names=TRUE, col.names = NA,
            x = iteratedRf.test.edgeR)
summaryRFedgeRTidy <- tidy(summary(iteratedRf.test.edgeR$testing.error.rate))
write.table(file = file.path(table.dirs[classifier], "RF_testing_for_edger_ordered_log2-normalised_variables_with_real_dataset.tsv"),
            quote=FALSE, sep="\t", row.names=TRUE, col.names = NA,
            x = summaryRFedgeRTidy)
pdf(file = file.path(dir.rFsubFigures, paste("RF_testing_with_edgeR_odered_with_log2-normalised_for_real_dataset", ".pdf", sep = "")),
    width =9 , height = 7 )
boxplot(iteratedRf.test.edgeR[,c("training.error.rate", "testing.error.rate")],
        ylim = c(0,1),
        main= "RF_testing_with_edgeR_odered_with_log2-normalised_for_real_dataset")
silence <- dev.off()

########################################################
## results after we applying DESeq2 DE and arranging the variables with the most significance to study
## the effect of such ordering onto the efficaincy of the RF classifier.

#### Results of RF real data set with DESseq2 ordaring variables
write.table(file = file.path(table.dirs[classifier], "Testing_RF_with_DESeq2_ordered_variables_with_log2-normalised_for_real_dataset.tsv"),
            quote=FALSE, sep="\t", row.names=TRUE, col.names = NA,
            x = iterated.test.DESeq2)
summaryRFDESeq2Tidy <- tidy(summary(iterated.test.DESeq2$testing.error.rate))
write.table(file = file.path(figure.dirs[classifier], "RF_summary_testing_with_DESeq2_ordered_variables_log2_normalised_for_real_dataset.tsv"),
            quote=FALSE, sep="\t", row.names=TRUE, col.names = NA,
            x = summaryRFDESeq2Tidy)
pdf(file = file.path(dir.figures, paste("RF_testing_with_DESeq2_ordered_variables_log2_normalised_for_real_dataset", ".pdf", sep = "")),
    width =9 , height = 5 )
boxplot(iterated.test.DESeq2[,c("training.error.rate", "testing.error.rate")],
        ylim = c(0,1),
        main= "Comparing training and testing error rate for DESeq2 ordaring variables with RF")
silence <- dev.off()


###########################################################################################################
###### Testing random Forest classifier performance with the DESeq2 and edgeR ordaring variables for #####
###### permuted data set                                                                             ######
###########################################################################################################
## we thus will using the ordered variables which are the most significance resulting from DEG,
## concering with padj.
if(parameters$compute){
  message.with.time("analysis RF classifier with DESeq2 and edgeR ordaring permuted data set , ", parameters$iterations, " iterations.")

  # DEG.DESeq2 <- DEGordering(loaded$countTable , loaded$classes, method = "DESeq2")
  # DEG.edgeR  <- DEGordering(loaded$countTable, loaded$classes, method = "edgeR")

  sorted.log2.transformed.edgeR <-log2norm$Counts[, DEG.DESeq2$geneOrder]
  sorted.log2.transformed.DESeq2 <- log2norm$Counts[, DEG.edgeR$geneOrder]

  iteratedRf.test.edgeRPermuted <- data.frame()
  iteratedRf.test.DESeq2Permuted<- data.frame()

  edgeR.error.table <- data.frame()
  DESeq2.error.table <- data.frame()


  for(i in 1:length(parameters$nb.variables)){
    varnb <- parameters$nb.variables[i]

    for(i in 1:parameters$iterations){
      test.edgeR <-MisclassificationEstimate(countTable = sorted.log2.transformed.edgeR[, 1:varnb], classes = sample(classes) ,classifier = "rf" )
      test.DESeq2 <-MisclassificationEstimate(countTable = sorted.log2.transformed.DESeq2[, 1:varnb], classes = sample(classes), classifier = "rf")
      iteratedRf.test.edgeRPermuted <- rbind(iteratedRf.test.edgeRPermuted, test.edgeR$stats)
      iteratedRf.test.DESeq2Permuted <- rbind(iteratedRf.test.DESeq2Permuted, test.DESeq2$stats)

    } # end for iteration
  } # end for number of variables
} else
{
  message.with.time("Skipping testing RF classifier with DESeq2 & edgeR ordered variables for log2-normalised permuted dataset")
}
## Hereafter we noticed that ordaring of variable has achieved much improvement in the calssification
## performing the Differential expression analysis between our two classes (Heparinised.blood & Bone marrow)
## will foster the efficiancy's classifier

## results after we applying edgeR DE and arranging the variables with the most significance to study
## the effect of such ordering onto the efficaincy of the RF classifier.
#iterated.test.edgeR

### Results of FR with permuted data set with edgeR ordaring variables
write.table(file = file.path(table.dirs[classifier], "Testing_RF_with_edgeR_ordered_with_log2-normalised_permuted_dataset.tsv"),
            quote=FALSE, sep="\t", row.names=TRUE, col.names = NA,
            x = iteratedRf.test.edgeRPermuted)
summaryRFedgeRTidy <- tidy(summary(iteratedRf.test.edgeRPermuted$testing.error.rate))
write.table(file = file.path(table.dirs[classifier], "Testing_RF_with_edgeR_ordered_with_log2-normalised_permuted_dataset.tsv"),
            quote=FALSE, sep="\t", row.names=TRUE, col.names = NA,
            x = summaryRFedgeRTidy)
pdf(file = file.path(figure.dirs[classifier], paste( "Testing_RF_with_edgeR_ordered_with_log2-normalised_permuted_dataset", ".pdf", sep = "")),
    width =9 , height = 7 )
boxplot(iteratedRf.test.edgeRPermuted[,c("training.error.rate", "testing.error.rate")],
        ylim = c(0,1),
        main= "Testing_RF_with_edgeR_ordered_with_log2-normalised_permuted_dataset")
silence <- dev.off()

########################################################
## results after we applying DESeq2 DE and arranging the variables with the most significance to study
## the effect of such ordering onto the efficaincy of the RF classifier.

#### Results of RF permuted data set with DESseq2 ordaring variables
write.table(file = file.path(table.dirs[classifier], "Testing_RF_with_DESeq2_ordered_with_log2-normalised_permuted_dataset.tsv"),
            quote=FALSE, sep="\t", row.names=TRUE, col.names = NA,
            x = iteratedRf.test.DESeq2Permuted)
summaryRFDESeq2Tidy <- tidy(summary(iteratedRf.test.DESeq2Permuted$testing.error.rate))
write.table(file = file.path(table.dirs[classifier], "Testing_RF_with_DESeq2_ordered_with_log2-normalised_permuted_dataset.tsv"),
            quote=FALSE, sep="\t", row.names=TRUE, col.names = NA,
            x = summaryRFDESeq2Tidy)
pdf(file = file.path(figure.dirs[classifier], paste("Testing_RF_with_DESeq2_ordered_with_log2-normalised_permuted_dataset", ".pdf", sep = "")),
    width =9 , height = 7 )
boxplot(iteratedRf.test.DESeq2Permuted[,c("training.error.rate", "testing.error.rate")],
        ylim = c(0,1),
        main= "Testing RF with DESeq2 ordered with log2-normalised permuted dataset")
silence <- dev.off()

######################################################################################
###### collecting the results
######################################################################################
testingErrorRateTable <- data.frame(
  "KNN.raw"          = KNNrawTestTable$testing.error.rate,
  "KNN.norm"         = KNNnormTestTable$testing.error.rate,
  "KNN.norm.log2"     = iteratedKnnTestingErrorRatesLog2Norm$testing.error.rate,
  "RF.raw"           = iteratedRfTestingErrorRates$testing.error.rate,
  "RF.norm"          = iteratedRfTestingErrorRatesNorm$testing.error.rate,
  "RF.norm.log2"          = iteratedRfTestingErrorRatesLog2Norm$testing.error.rate,
  "KNN.raw.perm"   = KNNrawTestTablePermLabels$testing.error.rate,
  "KNN.norm.perm" = iteratedKnnTestingErrorRatesLog2NormPermuted$testing.error.rate,
  "RF.raw.perm"   = iteratedRfTestingErrorRatesPermuted$testing.error.rate,

  "RF.norm.perm"  = iteratedRfTestingErrorRatesNormPermuted$testing.error.rate,
  "KNN.norm.log2.perm"  = iteratedKnnTestingErrorRatesLog2NormPermuted$testing.error.rate,
  "RF.norm.log2.perm"  = iteratedRfTestingErrorRatesLog2NormPermuted$testing.error.rate,

  # "Knn.Ordered.edgeR"  = iteratedKnn.test.edgeR[, "testing.error.rate"],
  # "Knn.Ordered.DESeq2"  = iteratedKnn.test.DESeq2[, "testing.error.rate"],
  # "Rf.Ordered.edgeR" =   iteratedRf.test.edgeR[, "testing.error.rate"],
  # "Rf.Ordered.DESeq2" =   iteratedRf.test.DESeq2[, "testing.error.rate"],
  #
  # "Knn.Ordered.edgeR.Permuted"  = iteratedKnn.test.edgeRPermuted[, "testing.error.rate"],
  # "Knn.Ordered.DESeq2.Permuted"  = iteratedKnn.test.DESeq2Permuted[, "testing.error.rate"],
  # "Rf.Ordered.edgeR.Permuted"   =   iteratedRf.test.edgeRPermuted[, "testing.error.rate"],
  # "Rf.Ordered.DESeq2.Permuted"   =   iteratedRf.test.DESeq2Permuted[, "testing.error.rate"]


  "RF.norm.perm"  = iteratedRfTestingErrorRatesNormPermuted$testing.error.rate
)
# View(testingErrorRateTable)

## MUSTAFA: TO DO

next.steps <- data.frame(
  "KNN.norm.log2.perm"  = iteratedKNNTestingErrorRatesLog2NormPermuted$testing.error.rate,
  "RF.norm.log2.perm"  = iteratedRFTestingErrorRatesLog2NormPermuted$testing.error.rate,
  KnnOrderededgeR  = iteratedKnn.test.edgeR[, "testing.error.rate"],
  KnnOrderedDESeq2  = iteratedKnn.test.DESeq2[, "testing.error.rate"],
  RfOrderededgeR =   iteratedRf.test.edgeR[, "testing.error.rate"],
  RfOrderedDESeq2 =   iteratedRf.test.DESeq2[, "testing.error.rate"],

  KnnOrderededgeRPermuted  = iteratedKnn.test.edgeRPermuted[, "testing.error.rate"],
  KnnOrderedDESeq2Permuted  = iteratedKnn.test.DESeq2Permuted[, "testing.error.rate"],
  RfOrderededgeRPermuted   =   iteratedRf.test.edgeRPermuted[, "testing.error.rate"],
  RfOrderedDESeq2Permuted   =   iteratedRf.test.DESeq2Permuted[, "testing.error.rate"]

)

write.table(file = file.path(
  dir.results,
  paste(sep="", "classifier_comparisons_testing_error_rate_R", parameters$iterations, ".tsv")),
  x = testingErrorRateTable,
  sep="\t", quote=FALSE, row.names = TRUE, col.names = NA)
#summary(testingErrorRateTable)




## Draw a box plot to compare the error rates for all the classifier tests
pdf(file = file.path(
  dir.figures,
  paste(sep="",
        parameters$recountID,
        "_error_rate_boxplot_R", parameters$iterations,
        ".pdf")), width = 8, height = 11)

save.margins <- par("mar") ## Save the original margin values
par(mar=c(3,10,4,1)) ## Set a wide left margin to display Y axis labels in the boxplot

testingErrorColors <- c(1:ncol(testingErrorRateTable))
names(testingErrorColors) <- names(testingErrorRateTable)
testingErrorColors[grep(names(testingErrorColors), pattern="KNN")] <- "#9900FF"
testingErrorColors[grep(names(testingErrorColors), pattern="RF")] <- "#00BB00"
testingErrorColors[grep(names(testingErrorColors), pattern="perm")] <- "grey"

boxplot(testingErrorRateTable, horizontal = TRUE, las=1, cex.axis=0.8, ylim=c(0,1),
        main=paste(sep="", "Misclassification rates\n", parameters$iterations, " iterations"),
        xlab="Misclassification rate",
        col=testingErrorColors)
abline(v=seq(0, 1, 0.05), col="#DDDDDD", lty="dashed")
abline(v=seq(0, 1, 0.1), col="#AAAAAA", lty="dashed")
par(mar=save.margins) ## Restore the normal margins for next drawings
silence<- dev.off()



## Draw a box plot to compare the error rates for all the classifier tests
pdf(file = file.path(
  dir.figures,
  paste(sep="",
        parameters$recountID,
        "_error_rate_boxplot_R", parameters$iterations,
        ".pdf")), width = 8, height = 11)

save.margins <- par("mar") ## Save the original margin values
par(mar=c(3,10,4,1)) ## Set a wide left margin to display Y axis labels in the boxplot

testingErrorColors <- c(1:ncol(testingErrorRateTable))
names(testingErrorColors) <- names(testingErrorRateTable)
testingErrorColors[grep(names(testingErrorColors), pattern="KNN")] <- "#9900FF"
testingErrorColors[grep(names(testingErrorColors), pattern="RF")] <- "#00BB00"
testingErrorColors[grep(names(testingErrorColors), pattern="perm")] <- "grey"

boxplot(testingErrorRateTable, horizontal = TRUE, las=1, cex.axis=0.8, ylim=c(0,1),
        main=paste(sep="", "Misclassification rates\n", parameters$iterations, " iterations"),
        xlab="Misclassification rate",
        col=testingErrorColors)
abline(v=seq(0, 1, 0.05), col="#DDDDDD", lty="dashed")
abline(v=seq(0, 1, 0.1), col="#AAAAAA", lty="dashed")
par(mar=save.margins) ## Restore the normal margins for next drawings
silence<- dev.off()
############################################################################################
### analysis of the impact of the number of ordered vaiables by DEG-p-vale into KNN classifier
### with log2-transformed normalised counts

testingKnnOrderedVariablesReal <- data.frame(
  NoOfVariales        = iteratedVarnb ,
  iteratedKnnTestedgeR = meanOfEachIteration["testing.error.rate"],
  iteratedKnnTestedgeRAll =    iteratedmeanOfEachIterationAll["testing.error.rate"]
)

write.table(file = file.path(dir., "impact_of_number_of_ordered_variables_into_KNN_all.tsv"),
            testingKnnOrderedVariablesAll, sep = "\t", quote = F, col.names = NA, row.names = T)


pdf(file = file.path(dir.figures, paste(parameters$recountID, "Impact_of_number_of_order_variables_into_KNN_all",parameters$iterations
                                        ,".pdf")), width = 9 , height = 7)

save.margins <- par("mar")
par(mar = c(3,15,4,1))
testingKnnOrderedColersAll <- c(1:ncol(testingKnnOrderedVariablesReal))
names(testingKnnOrderedColersAll)<-names(testingKnnOrderedVariablesReal)
testingKnnOrderedColersAll[grep(names(testingKnnOrderedColersAll), pattern = "Log2Norm")]<-"#9900FF"
testingKnnOrderedColersAll[grep(names(testingKnnOrderedColersAll), pattern = "allVariables")] <- "grey"

boxplot(testingKnnOrderedVariablesReal, horizontal = T, las=1, cex.axis=0.8, ylim= c(0,1),
        xlab= "Misclassification Rate",
        main= paste("Misclassification Rates\n",parameters$iterations, "iterations"),
        col= testingKnnOrderedColersAll)
abline(v= seq(0 , 1 , .05), col= "#DDDDDD", lty= "dashed")
abline(v=seq(0, 1 , .01),col= "#DDDDDD" ,lty= "dashed" )
par(mar = save.margins)
silence<- dev.off()
#################################################################################
#################### impact of each ordered variables
#################################################################################
testingKnnOrderedVariables <- data.frame(

  iteratedKnn.test.DESeq2.Log2Norm = iteratedKnn.test.DESeq2.Log2Norm[ c(1:10), "testing.error.rate"],
  iteratedKnn.test.DESeq2.Log2Norm = iteratedKnn.test.DESeq2.Log2Norm[c(11:20), "testing.error.rate"],
  iteratedKnn.test.DESeq2.Log2Norm = iteratedKnn.test.DESeq2.Log2Norm[c(21:30), "testing.error.rate"],
  iteratedKnn.test.DESeq2.Log2Norm = iteratedKnn.test.DESeq2.Log2Norm[c(31:40), "testing.error.rate"],
  iteratedKnn.test.DESeq2.Log2Norm = iteratedKnn.test.DESeq2.Log2Norm[c(41:50) , "testing.error.rate"],
  iteratedKnn.test.DESeq2.Log2Norm = iteratedKnn.test.DESeq2.Log2Norm[c(51:60), "testing.error.rate"],
  iteratedKnn.test.DESeq2.Log2Norm = iteratedKnn.test.DESeq2.Log2Norm[c(61:70), "testing.error.rate"],
  iteratedKnn.test.DESeq2.Log2Norm = iteratedKnn.test.DESeq2.Log2Norm[c(71:80), "testing.error.rate"],
  iteratedKnn.test.DESeq2.Log2Norm = iteratedKnn.test.DESeq2.Log2Norm[c(81:90), "testing.error.rate"],

  iteratedKnn.test.DESeq2.allVariables = iteratedKnn.test.DESeq2.allVariables[c(1:10), "testing.error.rate"],
  iteratedKnn.test.DESeq2.allVariables = iteratedKnn.test.DESeq2.allVariables[c(11:20), "testing.error.rate"],
  iteratedKnn.test.DESeq2.allVariables = iteratedKnn.test.DESeq2.allVariables[c(21:30), "testing.error.rate"],
  iteratedKnn.test.DESeq2.allVariables = iteratedKnn.test.DESeq2.allVariables[c(31:40), "testing.error.rate"],
  iteratedKnn.test.DESeq2.allVariables = iteratedKnn.test.DESeq2.allVariables[c(41:50), "testing.error.rate"],
  iteratedKnn.test.DESeq2.allVariables = iteratedKnn.test.DESeq2.allVariables[c(51:60), "testing.error.rate"],
  iteratedKnn.test.DESeq2.allVariables = iteratedKnn.test.DESeq2.allVariables[c(61:70), "testing.error.rate"],
  iteratedKnn.test.DESeq2.allVariables = iteratedKnn.test.DESeq2.allVariables[c(71:80), "testing.error.rate"],
  iteratedKnn.test.DESeq2.allVariables = iteratedKnn.test.DESeq2.allVariables[c(81:90), "testing.error.rate"]
)

write.table(file = file.path(dir.results, "impact_of_number_of_ordered_variables_into_KNN.tsv"),
            testingKnnOrderedVariablesAll, sep = "\t", quote = F, col.names = NA, row.names = T)


pdf(file = file.path(dir.figures, paste(parameters$recountID, "Impact_of_number_of_order_variables_into_KNN",parameters$iterations
                                        ,".pdf")), width = 9 , height = 7)

save.margins <- par("mar")
par(mar = c(3,15,4,1))
testingKnnOrderedColers <- c(1:ncol(testingKnnOrderedVariables))
names(testingKnnOrderedColers)<-names(testingKnnOrderedVariables)
testingKnnOrderedColers[grep(names(testingKnnOrderedColers), pattern = "Log2Norm")]<-"#9900FF"
testingKnnOrderedColers[grep(names(testingKnnOrderedColers), pattern = "allVariables")] <- "grey"

boxplot(testingKnnOrderedVariables, horizontal = T, las=1, cex.axis=0.8, ylim= c(0,1),
        xlab= "Misclassification Rate",
        main= paste("Misclassification Rates\n",parameters$iterations, "iterations"),
        col= testingKnnOrderedColers)
abline(v= seq(0 , 1 , .05), col= "#DDDDDD", lty= "dashed")
abline(v=seq(0, 1 , .01),col= "#DDDDDD" ,lty= "dashed" )
par(mar = save.margins)
silence<- dev.off()


