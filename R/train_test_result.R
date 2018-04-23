#' @title Class to handle the results of  Training-Testing evaluation for a given classifier
#' @author Jacques van Helden and Mustafa AbuElQumsan
#' @descripiton EachAn object of this class contains the primary results of a training-testing evaluation of a classifier
#' @param dataset the dataset on which the classifier ran
#' @param classifier the classifier used to produce this result
#' @param call  the function call used to produce this result
#'
#' @return an object of class TrainResult containing the following attributes
#' \itemize {
#'   \item ... MUSTAFA WILL DESCRIBE THE ATTRIBUTES ...
#' }

#'
#' @export
TrainTestResult <- function(dataset,
                            classifier,
                            #                           call,
                            iteration,
                            trainPredictedClasses,
                            testPredictedClasses) {


  if (!is(object = dataset, class2 = "DataTableWithTrainTestSets")) {
    stop("TrainTestResult() error. The dataset must belong to class DataTableWithTrainTestSets")
  }

  if (!is.integer(as.integer(iteration))
      || (iteration < 1)
      || (iteration > dataset$trainTestProperties$iterations)
  ){
    stop("TrainTestResult() error. The iteration must be a Natural number between 1 and the max number of iterations (", dataset$trainTestProperties$iterations, "). ")
  }

  ## Build a first version of the object based on passed parameters
  object <- structure(
    list(
      ID = dataset$ID,
      classLabels = dataset$classLabels,
      dataType = dataset$dataType,
      parameters = dataset$parameters,
      classifier = classifier,
      #      call = call, we
      iteration = iteration,
      trainIndex = dataset$trainTestProperties$trainIndices[[iteration]],
      testIndex = dataset$trainTestProperties$testIndices[[iteration]],
      trainPredictedClasses = trainPredictedClasses,
      testPredictedClasses = testPredictedClasses
    ),
    class="TrainTestResult")




  ################################################################################
  ## Compute testing misclassification reate
  object$test.contingency <- table(object$classLabels[object$testIndex], testPredictedClasses) ## Compute contingency table
  ## A misclassificaiton error is defined as an observation for which the predicted class differs from the known class
  object$testing.errors <- object$classLabels[object$testIndex] != object$testPredictedClasses
  object$testing.error.nb <- sum(object$testing.errors)
  object$testing.error.rate  <- object$testing.error.nb / length(object$testing.errors)


  ################################################################################
  ## Compute misclassification rate on the training set (=learning error)
  object$train.contingency <- table(object$classLabels[object$trainIndex], trainPredictedClasses) ## Compute contingency table
  ## A misclassificaiton error is defined as an observation for which the predicted class differs from the known class
  object$training.errors <- object$classLabels[object$trainIndex] != object$trainPredictedClasses
  object$training.error.nb <- sum(object$training.errors)
  object$training.error.rate  <- object$training.error.nb / length(object$training.errors)


  if (dataset$parameters$verbose) {
    message("\t\t\t", classifier,
            "; training error rate = ", signif(digits=3, object$training.error.rate ),
            "; testing error rate = ", signif(digits=3, object$testing.error.rate ))
  }
  # Gather single-value stats in a vector
  object$stats <- data.frame(n = iteration,
                             trainSize =  length(object$trainIndex),
                             trainingProportion = dataset$trainTestProperties$trainingProportion,

                             # testing.predicted.classes = testPredictedClasses,
                             testing.error.nb = object$testing.error.nb,
                             testing.error.rate = object$testing.error.rate,

                             #  training.predicted.classes = trainPredictedClasses,
                             training.error.nb= object$training.error.nb,
                             training.error.rate = object$training.error.rate)

  #cn <- as.vector(colnames(object$stats)) ## keep colnames for cbind

  ## Build a vector from the contingency table in order to return a object in the form ofa list that can easily be cas as vector
  # contingency.df <- as.data.frame.table(test.contingency)
  # cont.names <- paste(contingency.df[,2], contingency.df[,1], sep="_pred_")
  # object$stats <- cbind(object$stats, t(contingency.df[,3]))
  # colnames(object$stats) <- c(cn, cont.names )
  # object$test.contingency <- test.contingencytest.contingency

  object$trainingProportion <- dataset$trainTestProperties$trainingProportion
  object$trainSize <-length(object$trainIndex)
  object$testSize <-length(object$testIndex)

  return(object)
}



#' @export
heatmap.TrainTestResult <- function(object) {

  heatmap(object$test.contingency, Rowv = NA, Colv = NA,
          col = gray.colors(n = 256), scale  = "column",
          main = paste(sep="", object$classifier, " ", object$ID,
                                                                      "; it = ", object$iteration))

}



#' @export
print.TrainTestResult <- function(object) {
  cat("Train-test result\n")
  cat("\tTest misclassificationr rate = ", object$testing.error.rate, "\n")
  cat("\tTrain (learning) misclassificationr rate = ", object$training.error.rate, "\n")

}
