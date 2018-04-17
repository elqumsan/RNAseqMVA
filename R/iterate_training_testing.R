#' @title Iterate training/testing procedures with a given classifier and a given data type
#' @author Mustafa ABUELQUMSAN and Jacques van Helden
#' @description for sake of the accuracy and due to the error rate have computed from sampleing from the origin count data,
#' # it is better to compute the the error rate multiple time and then we would find the avarage for the error rate.
#' @param self an object of class DataTableWithTrainTestSets
#' @param classifier is the type of classifier that is used with repeated process.
#' @param permute is show if the class lable are permuted this for sake of the knowing the strength and weaknesses of the classifier
#' @param file.prefix in order to let us to save file from the IterateTrainingTesting
#'
#' @return
#' @import foreach
#' @import doParallel
#' @export
#'
#'

IterateTrainingTesting <- function (self, ...) {
  message("\tRunning IterateTrainingTesting() with object of class\t", paste( collapse  = ",",class(self) ) )
  UseMethod("IterateTrainingTesting", self)
 # return(self)
}


IterateTrainingTesting.DataTableWithClasses <- function (self, ...) {
  message("\tRunning IterateTrainingTesting() with object of class\t", "DataTableWithClasses")
  self <- NextMethod("IterateTrainingTesting",self)
}

IterateTrainingTesting.default <- function(self, ...){
  message("\tFinished from IterateTrainingTesting() with object of class\t", paste(collapse = ",", class(self)))
 # return(self)
}


################################################################
## Define a function to iterate over one classifier with one particular data type.
IterateTrainingTesting.DataTableWithTrainTestSets <- function (self,
                                                               classifier, # supported: knn or rf
                                                               permute = FALSE, # permute the class labels before running the test
                                                               file.prefix = NULL # prefix for the saved files. If not provided, will be automatically generated
) {


  startTime <- Sys.time()


  ## Check the class of the object
  if (!is(object = self, class2 = "DataTableWithTrainTestSets")) {
    stop("IterateTrainingTesting() only accepts objects of class DataTableWithTrainTestSets")
  }

  ## Get parameters from the passed object
  parameters <- self$parameters
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
  if (is.null(self$trainTestProperties$trainIndices)) {
    stop("IterateTrainingTesting()  train indices are not defined")
  } else {
    trainIndices <- self$trainTestProperties$trainIndices
  }


  ## Check that number of vectors in trainIndices corresponds to number of iterations
  if (length(trainIndices) != iterations) {
    stop("Invalid specification of trainIndices: should be a list of vectors, with the same number of vectors as the iterations (", iterations, ").")
  }


  message.with.time("\t", "IterateTrainingTesting()",
                    "\n\tID: ", self$ID,
                    "\n\tclassifier: ", classifier,
                    "\n\tdata type: ", self$dataType,
                    "\n\tvariable type: ", self$variablesType,
                    "\n\tTrain/test iterations: ",   parameters$iterations)


  ## Define file prefix is not specified in paramters
  if (is.null(file.prefix)) {
    file.prefix <- paste(sep="_", classifier, self$ID,  self$dataType, self$variablesType)
    file.prefix <- sub(pattern = " ", replacement = "_", x = file.prefix) ## Avoid spaces in file names
    if (permute) {
      file.prefix <- paste(sep="_", file.prefix, "permLabels")
    }
  }

  ## Define directory based on the method
  dir.create(parameters$dir$tablesDetail[classifier], recursive = TRUE, showWarnings = F)
  message("\tOutput directory for result tables: ", parameters$dir$tablesDetail[classifier])


  ## Iterate train/test cycles
  testTable <- data.frame() ## Instantiate the test table
  if (project.parameters$global$parallel) {
    message("\t", format(Sys.time(), "%Y-%m-%d_%H%M%S"), "\t",
            recountID, "\t", classifier,
            "\ttrain/test\t", project.parameters$global$iterations, " iterations with ", project.parameters$global$no_cores, " cores.")
    ## Run a foreach loop and get the result back in a data frame with rbind.
    testTable <- foreach(i = 1:iterations, .combine = rbind) %dopar%
      MisclassificationEstimate(self = self,
                                iteration = i,
                                classifier = classifier,
                                permute = permute)$stats
  } else {
    i <- 1 #iterations
    for (i in 1:iterations) {
      ## Permute class labels if required
      # computing the testing errors rate for the KNN classfier
      # trainIndex <- trainIndices[[i]]
      message("\t", format(Sys.time(), "%Y-%m-%d_%H%M%S"), "\t",
              recountID, "\t", classifier,
              "\ttrain/test iteration ", i , "/", iterations)


      oneTest <- MisclassificationEstimate(self = self,
                                           iteration = i,
                                           classifier = classifier,
                                           permute = permute)

      testTable <- rbind (testTable, oneTest$stats)
    }
  } ## end if project.parameters$parallel


  ## Save the result table for KNN training/testing evaluation
  testTable.file <- file.path(parameters$dir$tablesDetail[classifier], paste(sep="", file.prefix, ".tsv"))

  write.table(file=testTable.file,
              x = testTable,
              row.names = TRUE, col.names = NA, sep="\t", quote=FALSE)

  ## Export a summary about the error rate
  errorSummary.file <- file.path(parameters$dir$tablesDetail[classifier], paste(sep="", file.prefix, "_summary.tsv"))
  write.table(tidy(summary(testTable$testing.error.rate)),
              quote=FALSE, sep="\t", row.names=TRUE, col.names = NA,
              file = errorSummary.file)

  ## Plot a box plot of training versus testing error
  boxplot.file <- file.path(
    parameters$dir$figuresDetail[classifier],
    paste(sep="", file.prefix, "_R", iterations, "_learning_vs_test_error_boxplot.pdf"))

  pdf(file=boxplot.file, width = 7, height = 5)
  main <- paste(sep="", classifier, "; ", self$dataType, " counts; ", self$variablesType, " variables")
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
  elapsedTimeFile <- file.path(parameters$dir$tablesDetail[classifier], paste(sep="", file.prefix, "_elapsed_time.txt"))
  write(file = elapsedTimeFile, x = paste(startTime, endTime, elapsedTime))
  message("Elapsed Time file: ", elapsedTimeFile)
  NextMethod("IterateTrainingTesting", self)
  return(testTable)
}
