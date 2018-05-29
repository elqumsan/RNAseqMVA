#' @title Class to handle the results of Training-Testing evaluation for a given classifier
#' @author Jacques van Helden and Mustafa AbuElQumsan
#' @description  Each an object of this class contains the primary results of a training-testing evaluation of a classifier
#' @param dataset the dataset on which the classifier ran
#' @param classifier the classifier used to produce this result
#' @param iteration is the number of resampling for sake of much more precision to evaluate used classifier.
#' @param trainPredictedClasses are the classes that are predicted based on training the classifier with trainIndices.
#' @param testPredicatedClasses are the classes that are predicated based on training an classifier with testIndices.
#' @param call  the function call used to produce this result
#'
#' @return an object of class TrainTestResult containing the following attributes
#' \itemize {
#'   \item ID: is that recountID for the test dataset.
#'   \item classLabels: are labels (classes) for each sample in the targeted dataset.
#'   \item dataType: are the feature type (preprocessing created data) to analyize the imapct of pre-processing into efficiency of classifer, e.g. (filtered, norm, log2norm, log2norm-Pcs, ....)
#'   \item parameters: are all general parameter joind with each dataset from methods init_recount_id().
#'   \item classifier: is the classifier used to produce such results.
#'   \item iteration: is number of resampling to individuals in order to paly it with classifier.
#'   \item trainIndices: are the whole indeices for the trainset to be used as the trainset for the classifier ran.
#'   \item testIndices: are the whole indices for the testset to be paased to classifier as the tesset part.
#'   \item trainPredictedClasses: all the classes that are predicted based on training the classifier with trainset.
#'   \item testPredictedClasses: all the classes that are predicted based on training the classifier with testset.
#'   \item test.contingency: is (testing error) that is a table which contains in raw the actual classes on the testing set and in column predicted classes from classifier on the testing set.wherein that in diagolan is the correct classification and out of diagonal is miscalssifcation error rate for the testing set.
#'   \item testing.errors: are the class labels that are existed in the actuall classes for the testset and  which are not existed in the predicted classes from the trained classifier.
#'   \item testing.error.nb: is the sume of the testing.errors.
#'   \item testing.error.rate: is the training.error.nb divided by whole counts of the testing errors.
#'   \item train.contingency: is (learning error) the table contain in the raw actual class and in column predicted class and in diagonal is the corect classification and out of diagonal is misclassification errors.
#'   \item training.errors: are the class labels that are founded in the actuall classes for the trainset and which are not existed in the predicted classes from the trained classifier.
#'   \item training.error.nb: is the sum of the training.errors
#'   \item training.error.rate: is the training.error.nb divided by overall counts of the training errors.
#'   \item
#'   \item trainProportion:  the ratio of the trainset from all the dataset.
#'   \item trainSize:        is volume of the train size from all the dataset it is computed by multiple number of individuals in the trainProportion magnitude.
#'   \item testSize:         the test size "remained size" from all dataset after we remove the trainSize from targeted dataset
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
  ## Compute testing misclassification rate
  object$test.contingency <- table(object$classLabels[object$testIndex], testPredictedClasses) ## Compute contingency table
  ## A misclassificaiton error is defined as an observation for which the predicted class differs from the known class
  object$testing.errors <- object$classLabels[object$testIndex] != object$testPredictedClasses
  object$testing.error.nb <- sum(object$testing.errors)t
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
