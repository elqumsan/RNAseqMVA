################################################################
#' @title  RNA-Seq classifier evaluation to assesse the performance of the classifier
#' @author Mustafa ABUELQUMSAN and Jacques van Helden
#' @description  this script to evaluate and assesse the performance of the RNA-Seq classifier by
#' # Random sampling (random partitioning) estimation of the misclassification rate.
#'
#' @param dataset  an object of class DataTableWithTrainTestSets
#' @param iteration current iteration number (the MisclassificationEstimate function is typically called iteratively) and it called number of resampling.
#' @param classifier is a type of the classifier.
#' @param permute=FALSE permute the calss labels to measure misclassifciation rate without relevant learning
#' @return
#'  \itemize{
#'      \item dataset = it is data table one row for the individual (sample) and one culomn for the feature (gene).
#'      \item classifier = what is the tested classifier, KNN, RF, SVM and LDA.
#'      \item iteration = nomber of resampling.
#'      \item testPredictedClasses = these are the test classes which are predicted from classifier.
#'      \item trainPredictedClasses = these are the train classes whci are predicated from classifier.
#'     }
#'
#' @example
#' oneTest <- MisclassificationEstimate(dataTable, classLabels, trainingProportion = 2/3, classifier = "rf")
#'
#' @import doMC
#' @import class
#' @import e1071
#' @import MASS
#' @import stats
#' @import randomForest
#'
#' @export
MisclassificationEstimate <- function(dataset,
                                      iteration,
                                      classifier,
                                      permute=FALSE) {


  # require(doMC)
  # registerDoMC(cores = 5)

  # ##
  # if (is(dataset, class2 = "PCsWithTrainTestSets") ){
  #   # message("\t\tUsing PCs as features for object of class ", paste(collapse=", ", class(dataset)))
  #   dataTable <- t(dataset$x)
  # } else {
    # message("\t\tUsing dataTable as features for object of class ", paste(collapse=", ", class(dataset)))

  ## We transpose the dataTable
  transposedDataTable <- t(dataset$dataTable)
  # }
 #  iteration <- dataset$trainTestProperties$iterations

  ## Get parameters from the passed object
  parameters <- dataset$parameters

  ## Check **required** parameters
  for (p in c("verbose", "iterations")) {
    if (is.null(parameters[[p]])) {
      stop("Missing required parameter: '", p,
           "'.\n\tPlease check configuration file. ")
    } else {
      assign(p, parameters[[p]])
    }
  }


  ## AssignGet sample classes from the object
  classLabels <- dataset$classLabels
  if (permute) {
    classLabels <- sample(classLabels, replace = FALSE)
  }

  trainIndex <- dataset$trainTestProperties$trainIndices[[iteration]]
  testIndex <- dataset$trainTestProperties$testIndices[[iteration]]
  # intersect(trainIndex, testIndex)
  #
  if (sum(1:dataset$nbSamples != sort(union(trainIndex, testIndex))) > 0) {
    stop("Inconsistent training and testing indices: the union differs from Nb sample (",dataset$nbSamples,"). ")
  }


  if (classifier == "knn"){
    if (is.null(parameters$knn$k)) {
      stop("Missing required parameter: 'parameters$knn$k'.\n\tPlease check configuration file. ")
    } else {
      k <- parameters$knn$k
    }

    # require("class")

    ## Compute predicted classes on the test set
    testClassifier <- knn(train = transposedDataTable[trainIndex, ],
                          test = transposedDataTable[testIndex, ],
                          cl = classLabels[trainIndex],
                          k = k)
    testPredictedClasses <- as.vector(testClassifier)

    ## Compute predicted calsses on the training set (this will be used to compute learning error)
    trainPredictedClasses <- knn(train = transposedDataTable[trainIndex, ],
                                 test = transposedDataTable[trainIndex, ],
                                 cl = classLabels[trainIndex],
                                 k = k)


    ## TO DO: the following code thousl be written only once, after the list of if ... else if ... else if over classifiers.
    ## For each classifier, return two vectors
    ## testPredictedClasses and trainPredictedClasses
    ## All the rest is the same -> no code dreplication.


  } else if (classifier == "rf") {

    # require("randomForest")
    ## we need to tune our predictive model by using multiple workers "cores", such step to run our code through parallel
    ##  rather than sequentially technologies
    # require(randomForest)
    #require(doMC)
    #registerDoMC(cores = 5)

    ## Computing Testing errors for Random Forest

    ## Train the classifier with the training subset
    trainedClassifier  <- randomForest(
      x = transposedDataTable[trainIndex, ],
      y =  as.factor(classLabels[trainIndex]),
      xtest = transposedDataTable[testIndex,],  keep.forest = T)
    ## MUSTAFA: I think you don't use the xtest result after this,
    ## since you use predict() to predict the class of the testinf set.
    ## You should either keep xtest as randomForest parameter, but then
    ## get the tesing predicted  classes from the ranomForest() result,
    ## or remove this option here, since you use predict below.

    ## Predict class of the testing subset
    testPredictedClasses<- predict(trainedClassifier , transposedDataTable[testIndex,] )

    # ## Compute testing errors
    # test.contingency <- table(classLabels[testIndex], rf.test.prediction)
    # testing.errors <- classLabels[testIndex] != rf.test.prediction
    # testing.error.nb <- sum(testing.errors)
    # testing.error.rate <- testing.error.nb / length(testing.errors)

    ## Computing training errors (= learning errors)
    trainPredictedClasses <- predict(trainedClassifier , transposedDataTable[trainIndex,] )

    ## MUSTAFA, I (JvH) replaced the following command by(randomForest)
    ## by the previous one (predict) because you are running the training
    ## two times.
    ##
    ## Once the classifier is trained, it is simpler to reuse it for computing
    ## both the testing error rate and the training error rate.
    # rf.train.prediction  <- randomForest(x = transposedDataTable[trainIndex, ],
    #                                    y = as.factor( classLabels[trainIndex]),
    #                                    xtest = transposedDataTable[trainIndex,], importance = T, keep.forest = T)
    #
    # training.randsampling.fit <- predict(rf.train.prediction , transposedDataTable[trainIndex,] )


    ## Compute training errors (=learning errors)
    # training.contingency <- table(classLabels[trainIndex], rf.train.prediction )
    # training.errors <- classLabels[trainIndex] != rf.train.prediction
    # training.error.nb <-sum(training.errors)
    # training.error.rate <- training.error.nb / length(training.errors)
    #message(" contingency Table fo randomForest," ,  sep = " , " ,test.contingency[,3])


  } else if (classifier == "lda") {

    ## we need to tune our predictive model by using multiple workers "cores", such step to run our code through parallel
    ##  rather than sequentially technologies

#    registerDoMC(cores = 5)

    # require("MASS")

    ## Computing Testing errors for Linear Discriminant Analysis

    ## Train the classifier with the training subset
    trainedClassifier  <- lda(
      x = transposedDataTable[trainIndex, ],
      grouping = as.vector(classLabels[trainIndex]),
      CV=FALSE)

    ## Predict class of the testing subset
    testPredictedClasses <- predict(trainedClassifier , transposedDataTable[testIndex,])

    # test.contingency <- table(classLabels[testIndex], as.vector(lda.test.prediction$class))
    # testing.errors <- classLabels[testIndex] != lda.test.prediction$class
    # testing.error.nb <- sum(testing.errors)
    # testing.error.rate <- testing.error.nb / length(testing.errors)

    ## computing Training errors
    trainPredictedClasses <- predict(trainedClassifier , transposedDataTable[trainIndex,] )

    # training.contingency <- table(classLabels[trainIndex], as.vector(lda.train.prediction$class))
    # training.errors <- classLabels[trainIndex] != lda.train.prediction$class
    # training.error.nb <-sum(training.errors)
    # training.error.rate <- training.error.nb / length(training.errors)
    #message(" contingency Table fo randomForest," ,  sep = " , " ,test.contingency[,3])

    # result$trainingProportion <- trainingProportion
    # result$testing.error.nb <- testing.error.nb
    # result$testing.error.rate <-   testing.error.rate
    # result$training.error.nb <- training.error.nb
    # result$training.error.rate <-   training.error.rate



  } else if (classifier == "svm"){

  ## we need to tune our predictive model by using multiple workers "cores", such step to run our code through parallel
  ##  rather than sequentially technologies
    #registerDoMC(cores = 5)

    # require("e1071")

    ## Computing testing error for Suuport Vector Machine
    if (is.null(parameters$svm$kernel)) {
      parameters$svm$kernel <- "linear"
    }


    ## Train the classifier with the training suset
    trainedClassifier <- svm(x = transposedDataTable[trainIndex,],
                       y = as.factor(classLabels[trainIndex]),
                       type = parameters$svm$type,
                       scale = parameters$svm$scale,
                      # kernel = parameters$svm$kernel,
                       kernel = dataset$parameters$svm$kernel
                       # gamma = 1 ,
                       # cost = 10
    )

    ## predicting the classes of testing subset
    testPredictedClasses <- predict(trainedClassifier, transposedDataTable[testIndex,])
    length(testPredictedClasses)

    # test.contingency <- table(classLabels[testIndex], svm.test.prediction)
    # testing.errors <- classLabels[testIndex] != svm.test.prediction
    # testing.error.nb <- sum(testing.errors)
    # testing.error.rate <- testing.error.nb / length(testing.errors)

    ## Compute Training errors
    trainPredictedClasses <- predict(trainedClassifier, transposedDataTable[trainIndex,])
    length(trainPredictedClasses)

    # training.contingency <- table(classLabels[trainIndex], svm.train.prediction)
    # training.errors <- classLabels[trainIndex] != svm.train.prediction
    # training.error.nb <- sum(training.errors)
    # training.error.rate <- training.error.nb / length(training.errors)

  } else {
    stop(classifier, " is not a valid classifier. Supported: knn, rf, lda and svm" )

  }



  result <- TrainTestResult(dataset = dataset,
                            classifier = classifier,
                            iteration = iteration,
                            testPredictedClasses = testPredictedClasses,
                            trainPredictedClasses = trainPredictedClasses)

  return(result)
}

