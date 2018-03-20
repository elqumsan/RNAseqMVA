#################### Builde constractor for the selecting Training Sets #########
#' @title Export the object with train/test sets by stratification sampling
#' @author Mustafa ABUELQUMSAN and Jacques van Helden
#' @description sampling is done in each class separately in order to preserve the relative frequencies of classes in training and testing sets
#' @param self which much belong to the countTableWithClasses class.
#'
#' @export
  selectTrainingSets <- function(self, ...){
  message("\tExporting the object class  ", class(self), " object with train/test sets")
  UseMethod("selectTrainingSets", self)
}

## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ##
## TEMPORARILY HERE
## We define a method for stratified selection of training and testing sets,
## In a second time we will creat a new class names TrainTestCounttable, and this
## method will be attached to this class.

#' @title Random sampling of n training sets in a CountTableWithClasses.
#' @author Jacques van Helden and Mustafa AbuElQumsan
#' @description Select n random subsets for training. among the biological samples from a CountTableWithClasses.
#' @param self an object of the class CountTableWithClasses
#' @param stratified=TRUE if true, sampling is done in each class separately in order to preserve the relative frequencies of classes in training and testing sets.
#' @param iterations=parameters$iterations number of  train/test  iterations, which defines the number of independent sampled subsets
#' @param trainingProportion=parameters$trainingProportion proportion of samples to sample for each training set
#'
#' @export

  countTableWithTrainTestSets <- function( self,
                                           stratified=TRUE,
                                          iterations = parameters$iterations,
                                          trainingProportion = parameters$trainingProportion) {

  message.with.time("Selecting ", iterations, " training sets, with training proportion = ", trainingProportion)

  #### Check validity of the paraemters ####

  ## Check the class of input object
  if (!is(self, "countTableWithClasses")) {
    stop("selectStratifiedTrainingSets(): self parameter should belong to class countTableWithClasses. ")
  }
  ##  STRANGE: THIS RETURNS FALSE WHEREAS IT SHOULD B TRUE
  # isClass("countTableWithClasses")


  ## Trainng Proportion
  if ((trainingProportion < 0) || (trainingProportion > 1)) {
    stop("Training proportion must be a real number comprised between 0 and 1")
  }

  ## Instantiate the list with training indices
  trainIndices <- list()
  testIndices <- list()

  if (stratified) {
    ## Get class sizes
    trainSizePerClass <- round(self$samplesPerClass * trainingProportion)

    # testSizePerClass <- self$samplesPerClass - trainSizePerClass
    message("Stratified sampling among classes")
    print(as.data.frame(trainSizePerClass))
    # i <- 1
    for (i in 1:parameters$iterations) {
      trainIndices[[i]] <- vector()
      testIndices[[i]] <- vector()
      # c <- 1

      for (c in 1:self$nbClasses) {
        currentClass <- self$classNames[[c]]
        classSamples <- which (self$classLabels == currentClass)
        classTrain <- sample(classSamples, size = trainSizePerClass[[currentClass]], replace = FALSE)
        classTest <- setdiff(1:self$samplesPerClass[[c]],classTrain)
        trainIndices[[i]] <- append(trainIndices[[i]], classTrain)
        testIndices[[i]] <- append(testIndices[[i]], classTest)
        ## Check that the stratification  was correct
        ## table(self$classLabels[trainIndices[[i]]]) == trainSizePerClass

        #classTest <- setdiff(self$samplesPerClass[[c]], classTrain)
      }
    }
  } else {
    ## Sample the training sets irrespective of class membership
    n <- self$nbSamples
    trainSize <- round(trainingProportion * n)
    self$trainSize <- trainSize
    message("Class-independent sampling of training sets")
    for (i in 1:parameters$iterations) {
      trainIndices [[i]] <- sample(1:n, size = trainSize, replace = FALSE)
      testIndices [[i]] <- setdiff(1:n, trainIndices[[i]])
      #  View(as.data.frame.list(trainIndices))
    }
  }


  # ## Select testIndices as the complement of train indices
  # testIndices <- list()
  # for (i in 1:parameters$iterations) {
  #   ## MUTSAFA: DO IT
  #   setdiff()
  # }


  ## Add the attributes
  self$trainTestProperties <- list()
  self$trainTestProperties$iterations <- iterations
  self$trainTestProperties$trainingProportion <- trainingProportion
  self$trainTestProperties$trainSizePerClass <- trainSizePerClass
  self$trainTestProperties$trainIndices <- trainIndices
  self$trainTestProperties$testIndices <- testIndices
  self$trainTestProperties$stratified <- stratified

  class(self) <- unique(c( "countTableWithTrainTestSets", class(self)))

  return(self)
  message.with.time("Training set selection done")
}

summary.countTableWithTrainTestSets <- function(x){
  cat("\t\t countTableWithTrainTestSEts \n")
  cat("\tData Type       \t", x$dataType,"\n")
  cat("\tVariable Type      \t", x$variableType["all"], "\n")
  cat("\titeration          \t", x$trainTestProperties$iterations, "\n")
  cat("\ttraining Proportion   \t", x$trainTestProperties$trainingProportion, "\n")
  cat("\ttrain Size Per Class \t", x$trainTestProperties$trainSizePerClass, "\n")
  cat("\tstratified      \t", x$trainTestProperties$stratified, "\n")
  cat("\tclass     \t", class(x), "\n")

  #print()
}

print.countTableWithTrainTestSets <- function(x){
  summary.countTableWithTrainTestSets(x)
}
