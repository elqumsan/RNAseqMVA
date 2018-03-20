##### Iterate training/testing procedures with a given classifier and a given data type #####
#' @title Iterate training/testing procedures with a given classifier and a given data type
#' @author Mustafa ABUELQUMSAN and Jacques van Helden
#' @description for sake of the accuracy and due to the error rate have computed from sampleing from the origin count data,
#' # it is better to compute the the error rate multiple time and then we would find the avarage for the error rate.
#' @param countTable is the passed count data for the iterative process.
#' @param classes is the classes lable for the cout data.
#' @param classifier is the type of classifier that is used with repeated process.
#' @param data.type is data type that is used with iterative process, supported "raw","log2", "log2norm" and "log2norm.prcomp.centred"
#' @param iteration is the how much number it will be repeated the process.
#' @param variable.type this indicate for the number of variables that will used with iterative procedure.
#' @param trainingProportion is the ratio of trianing proportion.
#' @param trainIndices a list of vectors providing training sample indices that will be used at each iteration. The length of the list must correspond to the number of iteractions. This enables to compare different methods wilt using exactly the same training/testing sets between methods.
#' @param permute is show if the class lable are permuted this for sake of the knowing the strength and weaknesses of the classifier
#' @param file.prefix in order to let us to save file from the one experiment
#' @param verbose to write messages to indicate the progressing of tasks
#' @param k this parameter for the knn classifier to identify the number of neighbour in classification process.
#'
#' @return
#' @export
################################################################
## Define a function to iterate over one classifier with one particular data type.
one.experiment <- function (self,
                            classifier, # supported: knn or rf
                            iterations = parameters$iterations,
                            variable.type = "all", ## e.g;. "all", "top20", "top200", ...
                            trainingProportion, # ratio of training proportion
                          #  trainIndex, testIndex,
                            permute = FALSE, # permute the class labels before running the test
                            file.prefix = NULL, # prefix for the saved files. If not provided, will be automatically generated
                            verbose=TRUE,
                            k = parameters$knn$k ## For KNN only
) {


  startTime <- Sys.time()


  ## Check the class of the object
  if (!is(object = self, class2 = "countTableWithClasses")) {
    stop("one.experiment() only accepts objects of class countTableWithClasses")
  }

  ## Check the consistency between trainIndices and iterations
  if (!is.null(self$trainIndices)) {
      if (length(trainIndices) != iterations) {
        stop("Invalid specification of trainIndices: should be a list of vectors, with the same number of vectors as the iterations. ")
      }
  }

  if (!is(object = self, class2 = "countTableWithClasses")) {
    stop("one.experiment() only accepts objects of class countTableWithClasses")
  }

  message(format(Sys.time(), "%Y-%m-%d_%H%M%S"), "\t", classifier, " classifier (train vs test), ", self[["dataType"]], self[["varaible.type"]],  ",  ",variable.type, " variables, ",   parameters$iterations, " iterations.")
  testTable <- data.frame()
  # i <- 1



  ## Define file prefix is not specified in paramters
  if (is.null(file.prefix)) {
    file.prefix <- paste(sep="_", classifier, self$ID,  self$dataType, variable.type)
    if (permute) {
      file.prefix <- paste(sep="_", file.prefix, "permLabels")
    }
  }

  ## Define directory based on the method
  outDirectory <- table.dirs[classifier]
  dir.create(outDirectory, recursive = TRUE, showWarnings = F)
  message("Output directory: ", outDirectory)


  ## Iterate traint/test cycles
  for (i in 1:iterations) {
    ## Permute class labels if required
    if (permute) {
      classes <- sample(classes)
    }
    # trainIndex <- sample(trainIndex)
    # testIndex <- sample(testIndex)

    # computing the testing errors rate for the KNN classfier
    if (is.null(trainIndices)) {
      trainIndex <- NULL
    } else {
      trainIndex <- trainIndices[[i]]
    #  message("\t\ttrainIndex from trainIndices")
    }
    message("\t", format(Sys.time(), "%Y-%m-%d_%H%M%S"), "\t", classifier, " training/testing evaluation, iteration ", i , "/", iterations)

    oneTest <- MisclassificationEstimate(self,
                                          classifier = classifier,
                                          k=k,
                                          verbose=verbose)

    testTable <- rbind (testTable, oneTest$stats)
  }

  ## Save the result table for KNN training/testing evaluation
  testTable.file <- file.path(outDirectory, paste(sep="", file.prefix, ".tsv"))

  write.table(file=testTable.file,
              x = testTable,
              row.names = TRUE, col.names = NA, sep="\t", quote=FALSE)

  ## Export a summary about the error rate
  errorSummary.file <- file.path(outDirectory, paste(sep="", file.prefix, "_summary.tsv"))
  write.table(tidy(summary(testTable$testing.error.rate)),
              quote=FALSE, sep="\t", row.names=TRUE, col.names = NA,
              file = errorSummary.file)

  ## Plot a box plot of training versus testing error
  boxplot.file <- file.path(
    outDirectory,
    paste(sep="", file.prefix, "_R", iterations, "_learning_vs_test_error_boxplot.pdf"))

  pdf(file=boxplot.file, width = 7, height = 5)
  main <- paste(sep="", classifier, "; ", data.type, " counts; ", variable.type, " variables")
  if (permute) {
    main <- paste(sep="", main, "; permuted labels")
  }
  boxplot(testTable[, c("training.error.rate", "testing.error.rate")],
          ylim=c(0,1),
          main=main)
  silence <- dev.off()

    ## Compute elapsed time
  endTime <- Sys.time();
  elapsedTime <- endTime - startTime
  elapsedTimeFile <- file.path(outDirectory, paste(sep="", file.prefix, "_elapsed_time.txt"))
  write(file = elapsedTimeFile, x = paste(startTime, endTime, elapsedTime))
  message("Elapsed Time file: ", elapsedTimeFile)
  return(testTable)
}
