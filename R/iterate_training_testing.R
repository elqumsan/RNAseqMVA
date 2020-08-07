#' @title Iterate training/testing procedures with a given classifier and a given data type
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @description for sake of the accuracy and due to the error rate have computed from sampling from the origin count data,
#' it is better to compute the the error rate multiple time and then we would find the avarage for the error rate.
#' @param dataset an object of class DataTableWithTrainTestSets
#' @return an object which is Misclassification error rate for the specified number for resampling
#' \itemize{
#'      \item testTable: that is the table that is contains misclassification error rate for specified  number of resampling.
#'
#' }
#' @import foreach
#' @import doParallel
#' @export
#'
IterateTrainingTesting <- function(dataset) {
  message("\tRunning IterateTrainingTesting() with object of class\t", paste( collapse  = ",",class(dataset) ) )
  testTable <- UseMethod("IterateTrainingTesting", dataset)
  return(testTable)
}

#' @title Iterate training/testing procedures with a given classifier and a given data type for the object that belonge for DataTableWithClasses class
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @description for sake of the accuracy and due to the error rate have computed from sampleing from the origin count data,
#' it is better to compute the the error rate multiple time and then we would find the avarage for the error rate. with DataTableWithClasses's object.
#' @param dataset an object of class DataTableWithTrainTestSets
#' @return an object which is Misclassification error rate for the specified number for resampling
#'     \itemize{
#'         \item testTable: that is the table that is contains misclassification error rate for specified  number of resampling.
#'       }
#'
#' @export
IterateTrainingTesting.DataTableWithClasses <- function(dataset) {
  message("\tRunning IterateTrainingTesting() with object of class\t", "DataTableWithClasses")
  testTable <- NextMethod("IterateTrainingTesting", dataset)
  return(testTable)
}

#' @title printing report message to informaing so, Iterate training/testing procedures with a given classifier and a given data type for the object that belonge for DataTableWithClasses class has finished.
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @description printing report message to informaing so,for sake of the accuracy and due to the error rate have computed from sampleing from the origin count data,
#' it is better to compute the the error rate multiple time and then we would find the avarage for the error rate. with DataTableWithClasses's object.
#' @param dataset an object of class DataTableWithTrainTestSets
#'
#'
#' @export
IterateTrainingTesting.default <- function(dataset){
  message("\tFinished IterateTrainingTesting() with object of class\t", paste(collapse = ",", class(dataset)))
}


#' @title Iterate training/testing procedures with a given classifier and a given data type, with object has class DataTableWithTrainTestSets.
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @description for sake of the accuracy and due to the error rate have computed from sampleing from the origin count data,
#' it is better to compute the the error rate multiple time and then we would find the avarage for the error rate. with DataTableWithTrainTestSets's object
#' @param dataset an object of class DataTableWithTrainTestSets
#' @param classifier is the type of classifier that is used with repeated process.
#' @param permute is show if the class lable are permuted this for sake of the knowing the strength and weaknesses of the classifier
#' @param file.prefix prefix for files. If NULL, files will not be saved.
#' @param verbose=1 level of verbosity
#'
#' @return an object which is Misclassification error rate for the specified number for resampling
#'     \itemize{
#'         \item testTable: that is the table that is contains misclassification error rate for specified  number of resampling. for an object belonge to DataTableWithTrainTestSets.
#'       }
#' @export
IterateTrainingTesting.DataTableWithTrainTestSets <- function(
  dataset,
  classifier,
  permute = FALSE,
  file.prefix = NULL,
  verbose = 1
) {


  ## Check if file prefix has been correctly specified
  if (is.null(file.prefix)) {
    file.prefix <- paste(sep = "_", dataset$ID, dataset$dataType, classifier)
    if (permute) {
      file.prefix <- paste0(file.prefix, "_permLabels")
    }
  }
  startTime <- Sys.time()


  ## Check the class of the object
  if (!is(object = dataset, class2 = "DataTableWithTrainTestSets")) {
    stop("IterateTrainingTesting() only accepts objects of class DataTableWithTrainTestSets")
  }

  ## Get parameters from the passed dataset object
  parameters <- dataset$parameters

  ## Check required parameters
  for (p in c("verbose", "iterations")) {
    if (is.null(parameters[[p]])) {
      stop("Missing required parameter: '", p,
           "'.\n\tPlease check configuration file. ")
    } else {
      assign(p, parameters[[p]])
    }
  }


  ## Check the consistency between trainIndices and iterations
  if (is.null(dataset$trainTestProperties$trainIndices)) {
    stop("IterateTrainingTesting()  train indices are not defined")
  } else {
    trainIndices <- dataset$trainTestProperties$trainIndices
  }


  ## Check that number of vectors in trainIndices corresponds to number of iterations
  if (length(trainIndices) < iterations) {
    stop("Invalid specification of trainIndices: should be a list of vectors, at least as many vectors as the iterations (", iterations, ").")
  }

  if (verbose >= 1) {
    message.with.time("\t", "IterateTrainingTesting()",
                      "\n\tID: ", dataset$ID,
                      "\n\tfeature type: ", dataset$parameters$feature,
                      "\n\tdata type: ", dataset$dataType,
                      #                      "\n\tvariable type: ", dataset$variablesType,
                      "\n\tpermuted labels: ", permute,
                      "\n\tTrain/test iterations: ",   parameters$iterations,
                      "\n\tclassifier: ", classifier)

    if (classifier == "knn") {
      message("\tKNN k:  ", dataset$parameters$knn$k)
    } else if (classifier == "svm") {
      message("\tSVM kernel:  ", dataset$parameters$svm$kernel)
    }
  }


  ## Define directory based on the method
  # outParam <- outputParameters( dataset,classifier, permute)

  ## Iterate train/test cycles
  testTable <- data.frame() ## Instantiate the test table
  if (project.parameters$global$parallel) {
    if (verbose >= 1) {
      message("\t", format(Sys.time(), "%Y-%m-%d_%H%M%S"), "\t",
              recountID, "\t", classifier,
              "\ttrain/test\t", project.parameters$global$iterations, " iterations with ", project.parameters$global$no_cores, " cores.")
    }

    ## Run a foreach loop and get the result back in a data frame with rbind.
    testTable <- foreach(i = 1:parameters$iterations, .combine = rbind) %dopar%
      MisclassificationEstimate(dataset = dataset,
                                iteration = i,
                                classifier = classifier,
                                permute = permute)$stats

  } else {
    i <- 1 #iterations
    for (i in 1:parameters$iterations) {
      ## Permute class labels if required
      # computing the testing errors rate for the KNN classfier
      # trainIndex <- trainIndices[[i]]
      if (verbose >= 1) {
        message("\t", format(Sys.time(), "%Y-%m-%d_%H%M%S"), "\t",
                recountID, "\t", classifier,
                "\ttrain/test iteration ", i , "/", iterations)

      }

      oneTest <- MisclassificationEstimate(dataset = dataset,
                                           iteration = i,
                                           classifier = classifier,
                                           permute = permute)

      testTable <- rbind(testTable, oneTest$stats)
    }
  } ## end if project.parameters$parallel


  ## Save the result table for training/testing evaluation
  outDir <- parameters$dir$classifier_tables[classifier]
  dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
  testTable.file <- file.path(
    parameters$dir$classifier_tables[classifier],
    paste(sep = "", file.prefix, ".tsv"))

  message("\t\tSaving testing result table\t", testTable.file)
  write.table(file = testTable.file,
              x = testTable,
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

  ## Export a summary about the error rate
  errorSummary.file <- file.path(parameters$dir$classifier_tables[classifier], paste(sep = "", file.prefix, "_summary.tsv"))
  message("\t\tSaving testing error rates\t", errorSummary.file)
  write.table(summary(data.frame(testTable$testing.error.rate)),
              quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA,
              file = errorSummary.file)

  ## Plot a box plot of training versus testing error
  # message("\t\tBox plot")
  # boxplot.file <- file.path(
  #   parameters$dir$classifier_figures[classifier],
  #   paste(sep = "", file.prefix, "_R", iterations, "_learning_vs_test_error_boxplot.pdf"))
  #
  # pdf(file = boxplot.file, width = 7, height = 5)
  # main <- paste(sep = "", recountID, "\n", classifier, "; ", dataset$dataType, "; all variables")
  # if (permute) {
  #   main <- paste(sep = "", main, "; permuted labels")
  # }
  # boxplot(testTable[, c("training.error.rate", "testing.error.rate")],
  #         ylim = c(0,1),
  #         main = main)
  # silence <- dev.off(); rm(silence)

  ## Compute elapsed time
  endTime <- Sys.time();
  elapsedTime <- difftime(endTime, startTime, units = "secs")
  elapsedTimeFile <- file.path(parameters$dir$classifier_tables[classifier], paste(sep = "", file.prefix, "_elapsed_time.txt"))
  write(file = elapsedTimeFile, x = paste(sep = "\t", startTime, endTime, signif(digits = 3, elapsedTime)))
  message("\t\tElapsed Time file: ", elapsedTimeFile)
  NextMethod("IterateTrainingTesting", dataset)
  return(testTable)
}
