################################################################
#' @title  RNA-Seq classifier evaluation to assesse the performance of the classifier
#' @author Mustafa ABUELQUMSAN and Jacques van Helden
#' @description  this script to evaluate and assesse the performance of the RNA-Seq classifier by
#' # Random sampling (random partitioning) estimation of the misclassification rate.
#'
#' @param dataset  an object of class DataTableWithTrainTestSets
#' @iteration current iteration number (the MisclassificationEstimate function is typically called iteratively)
#' @param classifier is a type of the classifier
#' @param permute=FALSE permute the calss labels to measure misclassifciation rate without relevant learning
#' @example
#' oneTest <- MisclassificationEstimate(dataTable, classes, trainingProportion = 2/3, classifier = "rf")
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
  result <- list()

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
  ## Check required parameters
  for (p in c("verbose")) {
    if (is.null(parameters[[p]])) {
      stop("Missing required parameter: '", p,
           "'.\n\tPlease check configuration file. ")
    } else {
      assign(p, parameters[[p]])
    }
  }


  ## AssignGet sample classes from the object
  classes <- dataset$classLabels
  if (permute) {
    classes <- sample(classes, replace = FALSE)
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
    testPredictedClasses <- as.vector(knn(train = transposedDataTable[trainIndex, ],
                                          test = transposedDataTable[testIndex, ],
                                          cl = classes[trainIndex],
                                          k = k))

    ## Compute predicted calsses on the training set (this will be used to compute learning error)
    trainPredictedClasses <- knn(train = transposedDataTable[trainIndex, ],
                                 test = transposedDataTable[trainIndex, ],
                                 cl = classes[trainIndex],
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
    rf.trained  <- randomForest(
      x = transposedDataTable[trainIndex, ],
      y =  as.factor(classes[trainIndex]),
      xtest = transposedDataTable[testIndex,],  keep.forest = T)
    ## MUSTAFA: I think you don't use the xtest result after this,
    ## since you use predict() to predict the class of the testinf set.
    ## You should either keep xtest as randomForest parameter, but then
    ## get the tesing predicted  classes from the ranomForest() result,
    ## or remove this option here, since you use predict below.

    ## Predict class of the testing subset
    testPredictedClasses<- predict(rf.trained , transposedDataTable[testIndex,] )

    # ## Compute testing errors
    # test.contingency <- table(classes[testIndex], rf.test.prediction)
    # testing.errors <- classes[testIndex] != rf.test.prediction
    # testing.error.nb <- sum(testing.errors)
    # testing.error.rate <- testing.error.nb / length(testing.errors)

    ## Computing training errors (= learning errors)
    trainPredictedClasses <- predict(rf.trained , transposedDataTable[trainIndex,] )

    ## MUSTAFA, I (JvH) replaced the following command by(randomForest)
    ## by the previous one (predict) because you are running the training
    ## two times.
    ##
    ## Once the classifier is trained, it is simpler to reuse it for computing
    ## both the testing error rate and the training error rate.
    # rf.train.prediction  <- randomForest(x = transposedDataTable[trainIndex, ],
    #                                    y = as.factor( classes[trainIndex]),
    #                                    xtest = transposedDataTable[trainIndex,], importance = T, keep.forest = T)
    #
    # training.randsampling.fit <- predict(rf.train.prediction , transposedDataTable[trainIndex,] )


    ## Compute training errors (=learning errors)
    # training.contingency <- table(classes[trainIndex], rf.train.prediction )
    # training.errors <- classes[trainIndex] != rf.train.prediction
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
    lda.trained  <- lda(
      x = transposedDataTable[trainIndex, ],
      grouping = as.vector(classes[trainIndex]),
      CV=FALSE)

    ## Predict class of the testing subset
    testPredictedClasses <- predict(lda.trained , transposedDataTable[testIndex,])

    # test.contingency <- table(classes[testIndex], as.vector(lda.test.prediction$class))
    # testing.errors <- classes[testIndex] != lda.test.prediction$class
    # testing.error.nb <- sum(testing.errors)
    # testing.error.rate <- testing.error.nb / length(testing.errors)

    ## computing Training errors
    trainPredictedClasses <- predict(lda.trained , transposedDataTable[trainIndex,] )

    # training.contingency <- table(classes[trainIndex], as.vector(lda.train.prediction$class))
    # training.errors <- classes[trainIndex] != lda.train.prediction$class
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
    svm.trained <- svm(x = transposedDataTable[trainIndex,],
                       y = as.factor(classes[trainIndex]),
                       type = parameters$svm$type,
                       scale = parameters$svm$scale,
                      # kernel = parameters$svm$kernel,
                       kernel = dataset$parameters$svm$kernel
                       # gamma = 1 ,
                       # cost = 10
    )

    ## predicting the classes of testing subset
    testPredictedClasses <- predict(svm.trained, transposedDataTable[testIndex,])

    # test.contingency <- table(classes[testIndex], svm.test.prediction)
    # testing.errors <- classes[testIndex] != svm.test.prediction
    # testing.error.nb <- sum(testing.errors)
    # testing.error.rate <- testing.error.nb / length(testing.errors)

    ## Compute Training errors
    trainPredictedClasses <- predict(svm.trained, transposedDataTable[trainIndex,])

    # training.contingency <- table(classes[trainIndex], svm.train.prediction)
    # training.errors <- classes[trainIndex] != svm.train.prediction
    # training.error.nb <- sum(training.errors)
    # training.error.rate <- training.error.nb / length(training.errors)

  } else {
    stop(classifier, " is not a valid classifier. Supported: knn, rf, lda and svm" )

  }





  ################################################################################
  ## Compute testing misclassification reate
  test.contingency <- table(classes[testIndex], testPredictedClasses) ## Compute contingency table
  ## A misclassificaiton error is defined as an observation for which the predicted class differs from the known class
  testing.errors <- classes[testIndex] != testPredictedClasses
  testing.error.nb <- sum(testing.errors)
  testing.error.rate  <- testing.error.nb / length(testing.errors)


  ################################################################################
  ## Compute misclassification rate on the training set (=learning error)
  train.contingency <- table(classes[trainIndex], trainPredictedClasses) ## Compute contingency table
  ## A misclassificaiton error is defined as an observation for which the predicted class differs from the known class
  training.errors <- classes[trainIndex] != trainPredictedClasses
  training.error.nb <- sum(training.errors)
  training.error.rate  <- training.error.nb / length(training.errors)


  if (verbose) {
    message("\t\t\t", classifier,
            "; training error rate = ", signif(digits=3, training.error.rate ),
            "; testing error rate = ", signif(digits=3, testing.error.rate ))
  }
  # Gather single-value stats in a vector
  result$stats <- data.frame(n = iteration,
                             trainSize =  length(trainIndex),
                             trainingProportion = dataset$trainTestProperties$trainingProportion,

                           # testing.predicted.classes = testPredictedClasses,
                             testing.error.nb = testing.error.nb,
                             testing.error.rate = testing.error.rate,

                          #  training.predicted.classes = trainPredictedClasses,
                             training.error.nb= training.error.nb,
                             training.error.rate = training.error.rate)

  cn <- as.vector(colnames(result$stats)) ## keep colnames for cbind

  ## Build a vector from the contingency table in order to return a result in the form ofa list that can easily be cas as vector
  # contingency.df <- as.data.frame.table(test.contingency)
  # cont.names <- paste(contingency.df[,2], contingency.df[,1], sep="_pred_")
  # result$stats <- cbind(result$stats, t(contingency.df[,3]))
  # colnames(result$stats) <- c(cn, cont.names )
  # result$test.contingency <- test.contingencytest.contingency

  result$trainingProportion <- dataset$trainTestProperties$trainingProportion
  result$trainSize <-length(trainIndex)
  result$testing.error.nb <- testing.error.nb
  result$testing.error.rate <-   testing.error.rate

  return(result)
}

