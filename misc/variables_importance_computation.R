########### Computation the importace of all vaiables in raw count table and ordered it by the importance ########
countTable <-  rawCounts$Counts
classes <- loaded$classes

# n <- nrow(rawCounts$Counts) ## Number of observations (samples)
# trainingProportion <- 0.66
# train.size <- round(n * parameters$trainingProportion)
#
# ## Random selection of indices for the training set
# trainIndex <- sort(sample(1:n, size=train.size))
# ## Use remaining indices for the testing set
# testIndex <- setdiff(1:n, trainIndex)

rf.model  <- randomForest(
  x = countTable,
  y =  as.factor( classes),
  xtest = countTable, importance = T, keep.forest = T)

variable.importance <- importance(rf.model,type = 1,scale = F)
#variable.importance[,1]

ordered.varaible.importance <-order(variable.importance[,1],decreasing = T)

ordered.countTable.by.importace <-countTable[, ordered.varaible.importance]

sig.variables <- round(ncol(ordered.countTable.by.importace) * 0.75)
ordered.countTable.by.importance  <- ordered.countTable.by.importace[, 1:sig.variables]


classifier <- "svm" ## Default classifier for quick testing and debugging



if (parameters$identicalTrainTest) {
  ## New option: define all the train indices for all the iterations, in order to use the same training/testing sets between dfferent classifiers and data types
  trainIndices <- list()
  for(i in 1:parameters$iterations) {
    n <- nrow(rawCounts$Counts)
    train.size <- round(parameters$trainingProportion * n)
    trainIndices [[i]] <- sample(1:n, size = train.size, replace = FALSE)
  }
} else {
  ## First option: select different indices at each experiment
  trainIndices = NULL
}



# dim(counts)
# View(counts)

## Default for quick test without iterating over all cases
permute <- FALSE

if (parameters$compute) {

  for (classifier in parameters$classifiers) {

    ## List to store all results
    train.test.results.all.importance.varaibles <- list()
    #    train.test.results.all.PCs <- list()

    message.with.time("Train/test all computations with constant training proportion :",
                      signif(parameters$trainingProportion, digits = 3) )
    message.with.time("\tTrain/test, k=", parameters$knn$k, "; classifier=", classifier)



    #### Associate each analysis of real data with a permutation test ####
    for (permute in c(FALSE, TRUE)) {

      #### Run classifier with all variables (log2-transformed log counts) ####
      exp.prefix <-
        paste(sep = "_", classifier, parameters$recountID , parameters$data.type["V.importance"] , "allvars")
      if (permute) {
        exp.prefix <- paste(sep = "_", exp.prefix, perm.prefix)
      }# end if permuted class

      train.test.results.all.importance.varaibles[[exp.prefix]] <-
        one.experiment (
          countTable = as.data.frame(ordered.countTable.by.importance),
          classes = classes,
          trainIndices = trainIndices,
          # trainIndex = sample(log2norm$trainIndex),
          # testIndex = sample(log2norm$testIndex),
          data.type = parameters$data.types["V.importance"],
          classifier = classifier,
          #variable.type = variable.type,
          trainingProportion = parameters$trainingProportion,
          file.prefix = exp.prefix,
          permute = permute ,
          k = parameters$knn$k,
          verbose = parameters$verbose
        )



      #  } # end of iterative of Variable importance Numbers
    } # end for loop permutation


    #### Ploting the Misclassification Error rate for ordered varaibles by the importances of each variables #####  ####
    ErrorRateBoxPlot(experimentList = train.test.results.all.importance.varaibles,
                     classifier = classifier,
                     data.type = parameters$data.types["V.importance"],
                     main = paste(sep="",
                                  classifier, ": all importance variables,", "\n",
                                  parameters$recountID, ", ",
                                  parameters$iterations, " iterations, ","\n",
                                  data.type = "diverse data type"))


  } # end of loop over classifiers
} # end if statment computation

###############################################################################################
#### What is better to using ordered varaibles by the importances of each variables  ####
# ErrorRateBoxPlot(experimentList = train.test.results.all.importance.varaibles,
#                  classifier = classifier,
#                  data.type = parameters$data.types["V.importance"],
#                  main = paste(sep="",
#                               classifier, ": all importance variables,", "\n",
#                               parameters$recountID, ", ",
#                               parameters$iterations, " iterations, ","\n",
#                               data.type = "diverse data type"))


