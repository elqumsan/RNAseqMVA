#
# #studyCases$log2norm.prcomp.centred.scaled <- prcomp(t(na.omit(studyCases$log2norm$countTable)), center = TRUE, scale. = TRUE)
#
#
#
# message.with.time("finished running Principal Component analysis (PCA) for normalized log2 counts")
#
# ##### Visualization of the Principal Component analysis (PCA) to showcase impact of PCs #####
# message.with.time("Visualization of the Principal Component analysis (PCA) to showcase impact of PCs")
# source("misc/12_visualization_prcomp_analysis.R")


## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ##
## TEMPORARILY HERE
## We define a method for stratified selection of training and testing sets,
## In a second time we will creat a new class names TrainTestCounttable, and this
## method will be attached to this class.

#' @title Shorthand the real varaibles sets to the principal component.
#' @author Jacques van Helden and Mustafa AbuElQumsan
#' @description reduction all real varaibles in count Table to the n numbers from Principal component. among the biological samples from a CountTableWithClasses.
#' @param self an object of the class CountTableWithClasses

#' @export

PCsWithTrainTestSets <- function( self,
                                         stratified=TRUE,
                                         iterations = parameters$iterations,
                                         trainingProportion = parameters$trainingProportion) {

  message.with.time("Selecting ", iterations, " training sets, with training proportion = ", trainingProportion)

  #### Check validity of the paraemters ####

  ## Check the class of input object
  if (!is(self, "countTableWithTrainTestSets")  ) {
    stop("PCsWithTrainTestSets(): self parameter should belong to class countTableWithTrainTestSets. ")
  }
  ##  STRANGE: THIS RETURNS FALSE WHEREAS IT SHOULD B TRUE
  # isClass("countTableWithClasses")


  ## Trainng Proportion
  if ((trainingProportion < 0) || (trainingProportion > 1)) {
    stop("Training proportion must be a real number comprised between 0 and 1")
  }


  #### Compute principal components for normalized log2 counts ####
  message.with.time("Pre-processing by Principal Component analysis (PCA)")
  PCsProperties<- list()
  PCsProperties <- prcomp( t(na.omit(self$countTable)), center = TRUE, scale. = FALSE)
  PCsProperties$PCs <- PCsProperties$x

  ## Instantiate the list with training indices
  trainIndices <- list()
  testIndices <- list()

  #### iterate over permutation status ####
  pc.numbers <- c(2, 3, 4, 5, 6, 7,
                  seq(from=10, to=ncol(PCsProperties$PCs)-1, by = 10), ncol(PCsProperties$PCs))

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
      # trainIndices [[i]] <- sample(1:n, size = trainSize, replace = FALSE)
      # testIndices [[i]] <- setdiff(1:n, trainIndices[[i]])
      #
      trainSet <- sample(1:n, size = trainSize, replace = FALSE)
      testSet <- setdiff(1:n, trainSet)
      trainIndices[[i]] <-append(trainIndices[[i]], trainSet)
      testIndices[[i]] <- append(testIndices[[i]], testSet)

      #  View(as.data.frame.list(trainIndices))
    }
  }


  PCsProperties$ID <- self$ID
  PCsProperties$trainTestProperties <- self$trainTestProperties
  PCsProperties$sampleNames <- self$sampleNames
  PCsProperties$nbSamples<- self$nbSamples
  PCsProperties$variablesType <- "all"
  PCsProperties$dataType <-  "PCs_data"
  PCsProperties$classLabels <- self$classLabels
  PCsProperties$classNames <- self$classNames
  PCsProperties$nbClasses <- self$nbClasses
  PCsProperties$samplesPerClass <- self$samplesPerClass
  PCsProperties$sampleColors <- self$sampleColors
  PCsProperties$pc.numbers <- pc.numbers
  PCsProperties$variablesType <- "top_var"


  #self$trainTestProperties$stratified <- stratified

  #class(self) <- unique(c( "countTableWithTrainTestSets", class(self)))
  class(PCsProperties) <- unique(c( "PCsWithTrainTestSets", class(self)))

  return(PCsProperties)
  message.with.time("PCs Training set selection done")
}

summary.PCsWithTrainTestSets  <- function(x){
  cat("\t\t PCsWithTrainTestSets  \n")
  cat("\tData Type       \t", x$dataType,"\n")
  cat("\tVariables Type      \t", x$variablesType, "\n")
  cat("\titeration          \t", x$trainTestProperties$iterations, "\n")
  cat("\ttraining Proportion   \t", x$trainTestProperties$trainingProportion, "\n")
  cat("\ttrain Size Per Class \t", x$trainTestProperties$trainSizePerClass, "\n")
  cat("\tstratified      \t", x$trainTestProperties$stratified, "\n")
  cat("\tclass     \t", class(x), "\n")
  #print(x$trainTestProperties)
  cat("\n")

  #print()
}

print.countTableWithTrainTestSets <- function(x){
  summary.countTableWithTrainTestSets(x)
}
