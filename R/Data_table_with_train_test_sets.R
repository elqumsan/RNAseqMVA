#' @title Random sampling of n training sets in a DataTableWithClasses.
#' @author Jacques van Helden and Mustafa AbuElQumsan
#' @description Select n random subsets for training. among the biological samples from a DataTableWithClasses.
#' @param self an object of the class DataTableWithClasses
#' @param stratified=TRUE if true, sampling is done in each class separately in order to preserve the relative frequencies of classes in training and testing sets.
#' @param iterations=parameters$iterations number of  train/test  iterations, which defines the number of independent sampled subsets
#' @param trainingProportion=parameters$trainingProportion proportion of samples to sample for each training set
#'
#' @return object of class DataTableWithTrainTestSets containing the following attributes
#' \itemize{
#'  \item  iteration: which is number of resampling
#'   \item trainPropotion:      which is ratio of train volume from whole of dataset
#'   \item statified:           parameter which refer to stratification criteria would be included in analysis or not.
#'   \item trainSizePerClass:   which is specify the train volume for each class included in analysis.
#'   \item trainIndices:        which are whole indices for training process.
#'   \item testIndices:         which are whole indices for the testing process.
#'
#' }
#'
#' @export
DataTableWithTrainTestSets <- function(self,
                                        stratified=TRUE,
                                        iterations = parameters$iterations,
                                        trainingProportion = parameters$trainingProportion) {

  message.with.time("Building object of class DataTableWithTrainTestSets")

  #### Check validity of the paraemters ####

  ## Check the class of input object
  if (!is(self, "DataTableWithClasses")) {
    stop("selectStratifiedTrainingSets(): self parameter should belong to class DataTableWithClasses. ")
  }
  ##  STRANGE: THIS RETURNS FALSE WHEREAS IT SHOULD B TRUE
  # isClass("DataTableWithClasses")

  ## Assign test/train selection parameters as attributes of the object (to be passed to the attribute builder)
  self$trainTestProperties <- list()
  self$trainTestProperties$iterations <- iterations
  self$trainTestProperties$trainingProportion <- trainingProportion
  self$trainTestProperties$stratified <- stratified

  ## Add the class DataTableWithTrainTestSets
  class(self) <- unique(c(class(self), "DataTableWithTrainTestSets"))

  ## Build the attributes of the new object
  self <-buildAttributes(self)
#  self <- UseMethod( "buildAttributes", self)

  message("trainIndices length : ", length(self$trainTestProperties$trainIndices))
  message("\t\t returning from DataTableWithTrainTestSets()")
  return(self)
}


#' @title print a summary of an object belonging to class DataTableWithTrainTestSets
#' @author Mustafa AbuELQumsan and Jacques van Helden.
#' @param self is the targeted object for print its properities.
#' @description just print the summary of the object that is belonge to class DataTableWithTrainTestSets.
#' @return print the summary of such object that is belonge to DataTableWithTrainTestSets class by utilizing generic function.
#' @export
summary.DataTableWithTrainTestSets <- function(x){
  cat("\t\t DataTableWithTrainTestSEts \n")
  cat("\tData Type       \t", x$dataType,"\n")
  cat("\tVariables Type      \t", x$variablesType, "\n")
  cat("\titeration          \t", x$trainTestProperties$iterations, "\n")
  cat("\ttraining Proportion   \t", x$trainTestProperties$trainingProportion, "\n")
  cat("\ttrain Size Per Class \t", x$trainTestProperties$trainSizePerClass, "\n")
  cat("\tstratified      \t", x$trainTestProperties$stratified, "\n")
  cat("\tclass     \t", class(x), "\n")
  #print(x$trainTestProperties)
  cat("\n")

  NextMethod("summary", self)
  #print()
}

#' @title print of an object belonging to class DataTableWithTrainTestSets.
#' @author Mustafa AbuELQumsan and Jacques van Helden.
#' @param self is the targeted object for print its properities.
#' @description just print the summary of the object that is belonge to class DataTableWithTrainTestSets.
#' @return print the summary of such object that is belonge to DataTableWithTrainTestSets class by utilizing generic function.
#' @export
print.DataTableWithTrainTestSets <- function(x) {
  summary(x)
}




#' @title build the attributes that depend on the count table
#' @author Mustafa AbuELQumsan and Jacques van Helden
#' @description build special attributes for an object wherein to be initialized to belonge to DataTableWithTrainTestSets
#' @param self: is an object which we stive to initialize it to be has membership of DataTableWithTrainTestSets.
#' @return  an object that is became has the attributes of the class DataTableWithTrainTestSets
#'  \itemize{
#'  \item  iteration: which is number of resampling
#'   \item trainPropotion:      which is ratio of train volume from whole of dataset
#'   \item statified:           parameter which refer to stratification criteria would be included in analysis or not.
#'   \item trainSizePerClass:   which is specify the train volume for each class included in analysis.
#'   \item trainIndices:        which are whole indices for training process.
#'   \item testIndices:         which are whole indices for the testing process.
#'
#' }
#'
#' @export
buildAttributes.DataTableWithTrainTestSets <- function(self) {
  message("\tBuilding class-specific attributes for DataTableWithTrainTestSets\t", self$ID)

  ## Check stratified
  stratified <- self$trainTestProperties$stratified
  if (is.null(stratified)) {
    stop("Cannot build attributes for DataTableWithTrainTestSets because the attribute 'stratified' is null (should be Boolean)")
  }

  ## Check iterations
  iterations <- self$trainTestProperties$iterations
  if (is.null(iterations)) {
    stop("Cannot build attributes for DataTableWithTrainTestSets because the attribute 'iterations' is null (should be a strictly positive Natural)")
  }
  if (!is.integer(iterations) || (iterations < 1)) {
    stop("Invalid value for iterations (", iterations, "). Must be a strictly positive Natural number. ")

  }

  ## Check training proportion
  trainingProportion <- self$trainTestProperties$trainingProportion
  if (is.null(trainingProportion)) {
    stop("Cannot build attributes for DataTableWithTrainTestSets because the attribute 'trainingProportion' is null (should be a Real number)")
  }
  if (!is.numeric(trainingProportion) || (trainingProportion <= 0) || (trainingProportion >= 1)) {
    stop("Invalid value for trainingProportion (", trainingProportion, "). Must be a Real number between 0 and 1. ")

  }


  ## Instantiate the list with training indices
  message("\t\tSelecting ", iterations, " training sets, with training proportion = ", trainingProportion)
  trainIndices <- list()
  testIndices <- list()

  if (stratified) {
    ## Get class sizes
    trainSizePerClass <- round(self$samplesPerClass * trainingProportion)

    # testSizePerClass <- self$samplesPerClass - trainSizePerClass
    message("\t\tStratified sampling among classes")
    # print(as.data.frame(trainSizePerClass))
    # i <- 1
    for (i in 1:parameters$iterations) {
      trainIndices[[i]] <- vector()
      testIndices[[i]] <- vector()
      # c <- 1
      for (c in 1:self$nbClasses) {
        currentClass <- self$classNames[[c]]
        classSamples <- which (self$classLabels == currentClass)
        classTrain <- sample(classSamples, size = trainSizePerClass[[currentClass]], replace = FALSE)
        trainIndices[[i]] <- as.vector(append(trainIndices[[i]], classTrain))

        classTest <- setdiff(classSamples, classTrain)
        testIndices[[i]] <- append(testIndices[[i]], classTest)
        ## Check that the stratification  was correct
        ## table(self$classLabels[trainIndices[[i]]]) == trainSizePerClass

        #classTest <- setdiff(self$samplesPerClass[[c]], classTrain)
      }
      # length(trainIndices[[i]])
      # length(testIndices[[i]])
    }
  } else {
    ## Sample the training sets irrespective of class membership
    n <- self$nbSamples
    trainSize <- round(trainingProportion * n)
    self$trainSize <- trainSize
    message("\t\tClass-independent sampling of training sets")
    # i <- 1
    for (i in 1:parameters$iterations) {
      trainIndices[[i]] <- vector()
      testIndices[[i]] <- vector()
      trainSet <- sample(1:n, size = trainSize, replace = FALSE)
      testSet <- setdiff(1:n, trainSet)
      trainIndices[[i]] <-append(trainIndices[[i]], trainSet)
      testIndices[[i]] <- append(testIndices[[i]], testSet)
      #  View(as.data.frame.list(trainIndices))
    }
  }

  ## Sort train ant test indices (for the sake of readability)
  for (i in 1:parameters$iterations) {
    trainIndices[[i]] <- sort (trainIndices[[i]])
    testIndices[[i]] <- sort (testIndices[[i]])
  }

  ## Check that the sum of lengths for testindices and trainindices sum up to the length of the dataset
  if (sum(unlist(lapply(trainIndices, length)) + unlist(lapply(testIndices, length)) != self$nbSamples) > 0) {
    stop("Error with DataTableWithTrainTestSets(): incorrect lengths of trainIndices and testIndices ")
  }
  message.with.time("\tTraining set selection done")

  ## Convert trainIndices and testIndices from lists to
  ## data frames, in order to facilitate their further usage.
  # trainIndices <- as.data.frame.list(trainIndices)
  # colnames(trainIndices) <- 1:ncol(trainIndices)
  # rownames(trainIndices) <- 1:nrow(trainIndices)
  # testIndices <- as.data.frame.list(testIndices)
  # colnames(testIndicesgf) <- 1:ncol(testIndices)
  # rownames(testIndices) <- 1:nrow(testIndices)

  # ## Select testIndices as the complement of train indices
  # testIndices <- list()
  # for (i in 1:parameters$iterations) {
  #   ## MUTSAFA: DO IT
  #   setdiff()
  # }


  ## Add the attributes
  self$trainTestProperties$trainSizePerClass <- trainSizePerClass
  self$trainTestProperties$trainIndices <- trainIndices
  self$trainTestProperties$testIndices <- testIndices

#  print (self$trainTestProperties)
  self <- NextMethod("buildAttributes", self)


#  message("trainIndices length : ", length(self$trainTestProperties$trainIndices))
#  message("\t\treturning from buildAttributes.DataTableWithTrainTestSets()")
  return(self)
}
