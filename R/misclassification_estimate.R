################################################################
#'@title  RNA-Seq classifier evaluation to assesse the performance of the classifier
#'@author Jacques van Helden and Mustafa ABUELQUMSAN
#'@description  this script to evaluate and assesse the performance of the RNA-Seq calssifier by
#'# Random sampling (random partitioning) estimation of the misclassification rate.
#'
#'@param countTable  this data frame for RNA-Seq data which contains one row for indiviual and one column for variable
#'@param classes  such is vector for our case for classes,
#'@param k is number of neighbours passed to classifier
#'@param trainingproportion is the ratio of the training suset from the whole data set
#'@param calssfier is specie of the classifier
#'@example
#' oneTest <- MisclassificationEstimate(countTable, classes, trainingProportion = 2/3, classifier = "rf")
#'
#'

MisclassificationEstimate <- function(countTable, classes,
                                      trainingProportion = 2/3,
                                      classifier = "knn",
                                      verbose = FALSE,
                                      k= 3) {
  result <- list()
  n <- nrow(countTable) ## Number of observations (samples)
  train.size <- round(n * trainingProportion)

  ## Random selection of indices for the training set
  trainIndex <- sort(sample(1:n, size=train.size))
  ## Use remaining indices for the testing set
  testIndex <- setdiff(1:n, trainIndex)

  if (classifier == "knn"){

    ## we need to tune our predictive model by using multiple workers "cores", such step to run our code through parallel
    ##  rather than sequentially technologies
        library(doMC)
        registerDoMC(cores = 5)
    ## Compute testing errors
    randsampling.fit <- knn(train = countTable[trainIndex, ],
                            test = countTable[testIndex, ],
                            cl = classes[trainIndex],
                            k = k)

    test.contingency <- table(classes[testIndex], randsampling.fit) ## Compute contingency table
    ## A misclassificaiton error is defined as an observation for which the predicted class differs from the known class
    testing.errors <- classes[testIndex] != randsampling.fit
    testing.error.nb <- sum(testing.errors)

    testing.error.rate  <- testing.error.nb / length(testing.errors)
    #    message("knn, k=", k, "; testing error = ", signif(digits=3, testing.error.rate ))

    ## Compute training error ( = learning error)
    train.fit <- knn(train = countTable[trainIndex, ],
                     test = countTable[trainIndex, ],
                     cl = classes[trainIndex],
                     k = k)

    train.contingency <- table(classes[trainIndex], train.fit) ## Compute contingency table
    ## A misclassificaiton error is defined as an observation for which the predicted class differs from the known class
    training.errors <- classes[trainIndex] != train.fit
    training.error.nb <- sum(training.errors)

    training.error.rate  <- training.error.nb / length(training.errors)

  } else if (classifier == "rf") {

    ## we need to tune our predictive model by using multiple workers "cores", such step to run our code through parallel
    ##  rather than sequentially technologies

    library(doMC)
    registerDoMC(cores = 5)

    ## Computing Testing errors for Random Forest

    ## Train the classifier with the training subset
    rf.trained  <- randomForest(
      x = countTable[trainIndex, ],
      y = as.factor( classes[trainIndex]),
      xtest = countTable[testIndex,], importance = T, keep.forest = T)
    ## MUSTAFA: I think you don't use the xtest result after this,
    ## since you use predict() to predict the class of the testinf set.
    ## You should either keep xtest as randomForest parameter, but then
    ## get the tesing predicted  classes from the ranomForest() result,
    ## or remove this option here, since you use predict below.

    ## Predict class of the testing subset
    rf.test.prediction <- predict(rf.trained , countTable[testIndex,] )

    ## Compute testing errors
    test.contingency <- table(classes[testIndex], rf.test.prediction)
    testing.errors <- classes[testIndex] != rf.test.prediction
    testing.error.nb <- sum(testing.errors)
    testing.error.rate <- testing.error.nb / length(testing.errors)

    ## Computing training errors (= learning errors)
    rf.train.prediction <- predict(rf.trained , countTable[trainIndex,] )

    ## MUSTAFA, I (JvH) replaced the following command by(randomForest)
    ## by the previous one (predict) because you are running the training
    ## two times.
    ##
    ## Once the classifier is trained, it is simpler to reuse it for computing
    ## both the testing error rate and the training error rate.
    # rf.train.prediction  <- randomForest(x = countTable[trainIndex, ],
    #                                    y = as.factor( classes[trainIndex]),
    #                                    xtest = countTable[trainIndex,], importance = T, keep.forest = T)
    #
    # training.randsampling.fit <- predict(rf.train.prediction , countTable[trainIndex,] )


    ## Compute training errors (=learning errors)
    training.contingency <- table(classes[trainIndex], rf.train.prediction )
    training.errors <- classes[trainIndex] != rf.train.prediction
    training.error.nb <-sum(training.errors)
    training.error.rate <- training.error.nb / length(training.errors)
    #message(" contingency Table fo randomForest," ,  sep = " , " ,test.contingency[,3])


  } else if (classifier == "lda") {

    ## we need to tune our predictive model by using multiple workers "cores", such step to run our code through parallel
    ##  rather than sequentially technologies

    library(doMC)
    registerDoMC(cores = 5)

    library("MASS")

    ## Computing Testing errors for Linear Discriminant Analysis

    ## Train the classifier with the training subset
    lda.trained  <- lda(
      x = countTable[trainIndex, ],
      grouping = as.vector(classes[trainIndex]),
      CV=FALSE)

    ## Predict class of the testing subset
    lda.test.prediction <- predict(lda.trained , countTable[testIndex,])

    test.contingency <- table(classes[testIndex], as.vector(lda.test.prediction$class))
    testing.errors <- classes[testIndex] != lda.test.prediction$class
    testing.error.nb <- sum(testing.errors)
    testing.error.rate <- testing.error.nb / length(testing.errors)

    ## computing Training errors
    lda.train.prediction <- predict(lda.trained , countTable[trainIndex,] )
    training.contingency <- table(classes[trainIndex], as.vector(lda.train.prediction$class))
    training.errors <- classes[trainIndex] != lda.train.prediction$class
    training.error.nb <-sum(training.errors)
    training.error.rate <- training.error.nb / length(training.errors)
    #message(" contingency Table fo randomForest," ,  sep = " , " ,test.contingency[,3])

    # result$trainingProportion <- trainingProportion
    # result$testing.error.nb <- testing.error.nb
    # result$testing.error.rate <-   testing.error.rate
    # result$training.error.nb <- training.error.nb
    # result$training.error.rate <-   training.error.rate



  } else if (classifier == "svm"){

  ## we need to tune our predictive model by using multiple workers "cores", such step to run our code through parallel
  ##  rather than sequentially technologies
    library(doMC)
    registerDoMC(cores = 5)

    require("e1071")

    ## Computing testing error for Suuport Vector Machine
    if (is.null(parameters$svm$kernel)) {
      parameters$svm$kernel <- "radial"
    }


    ## Train the classifier with the training suset
    svm.trained <- svm(x = countTable[trainIndex,] ,
                       y = as.factor(classes[trainIndex]),
                       type = parameters$svm$type,
                       scale = parameters$svm$scale,
                       kernel = parameters$svm$kernel
                       # gamma = 1 ,
                       # cost = 10
    )

    ## predicting the classes of testing subset
    svm.test.prediction <- predict(svm.trained, countTable[testIndex,])

    test.contingency <- table(classes[testIndex], svm.test.prediction)
    testing.errors <- classes[testIndex] != svm.test.prediction
    testing.error.nb <- sum(testing.errors)
    testing.error.rate <- testing.error.nb / length(testing.errors)

    ## Compute Training errors
    svm.train.prediction <- predict(svm.trained, countTable[trainIndex,])

    training.contingency <- table(classes[trainIndex], svm.train.prediction)
    training.errors <- classes[trainIndex] != svm.train.prediction
    training.error.nb <- sum(training.errors)
    training.error.rate <- training.error.nb / length(training.errors)

  } else {
    stop(classifier, " is not a valid classifier. Supported: knn, rf, lda and svm" )

  }

  if (verbose) {
    message(classifier,
            "; training error rate = ", signif(digits=3, training.error.rate ),
            "; testing error rate = ", signif(digits=3, testing.error.rate ))
  }


  ## Gather single-value stats in a vector
  result$stats <- data.frame(n = n,
                             train.size = train.size,
                             trainingProportion = trainingProportion,
                             testing.error.nb = testing.error.nb,
                             testing.error.rate = testing.error.rate,
                             training.error.nb= training.error.nb,
                             training.error.rate = training.error.rate)
  cn <- as.vector(colnames(result$stats)) ## keep colnames for cbind

  ## Build a vector from the contingency table in order to return a result in the form ofa list that can easily be cas as vector
  # contingency.df <- as.data.frame.table(test.contingency)
  # cont.names <- paste(contingency.df[,2], contingency.df[,1], sep="_pred_")
  # result$stats <- cbind(result$stats, t(contingency.df[,3]))
  # colnames(result$stats) <- c(cn, cont.names )
  # result$test.contingency <- test.contingencytest.contingency

  result$trainingProportion <- trainingProportion
  result$testing.error.nb <- testing.error.nb
  result$testing.error.rate <-   testing.error.rate

  return(result)
}

