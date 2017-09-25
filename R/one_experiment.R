
################################################################
## Define a function to iterate over one classifier with one particular data type. 
one.experiment <- function (countTable, # count table (may be normalized or not, log2 or not)
                            classes, ## Class labels attached to each sample
                            classifier, # supported: knn or rf
                            data.type, # e.g. "raw", "normP75", "log2normP75"
                            iterations = parameters$iterations,
                            variable.type = "all", ## e.g;. "all", "top20", "top200", ...
                            trainingProportion, # ratio of training proportion
                            permute = FALSE, # permute the class labels before running the test
                            file.prefix = NULL, # prefix for the saved files. If not provided, will be automatically generated
                            verbose=FALSE,
                            k=3 ## For KNN only
) {
  
  
  startTime <- Sys.time();
  
  message(format(Sys.time(), "%Y-%m-%d_%H%M%S"), "\t", classifier, " classifier (train vs test), ", data.type, " counts, ", variable.type, " variables, ", parameters$iterations, " iterations.")
  testTable <- data.frame()
  # i <- 1
  
  
  ## Define file prefix is not specified in paramters
  if (is.null(file.prefix)) {
    file.prefix <- paste(sep="_", classifier, data.type, variable.type)
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
    
    # computing the testing errors rate for the KNN classfier
    message("\t", format(Sys.time(), "%Y-%m-%d_%H%M%S"), "\t", classifier, " training/testing evaluation, iteration ", i , "/", iterations)
    oneTest <- MisclassificationEstimate(countTable, classes , trainingProportion , classifier = classifier, k=k, verbose=verbose)
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
