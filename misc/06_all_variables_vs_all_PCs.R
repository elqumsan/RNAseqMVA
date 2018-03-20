
##### All variables versus all PCs. #####
##
## QUESTION: is it better to use the PCAs-transformed data, and, if so, is it better to use a subset of the first components or all the components ?
## For the time being we test this with only one classifier (KNN, default k)  but we will come back to it with other classifiers later.
## IMPORTANT NOTE : i would like to pay your intention for here we should take " data.type, so that we will not give the user to
## choose the data.type in return we will pass the data.type for each experiment.
## Choice of the classifier

# classifier <- "svm" ## Default classifier for quick testing and debugging
## Choice of the Counts
# data.type <- "log2norm.prcomp.centred"
# data.type <- "log2norm"

#
# if (parameters$identicalTrainTest) {
#   ## Define all the train indices for all the iterations, in order to use the
#   ## same training/testing sets between dfferent classifiers and data types
#   trainIndices <- list()
#   n <- nrow(log2norm$Counts)
#   train.size <- round(parameters$trainingProportion * n)
#   for(i in 1:parameters$iterations) {
#     trainIndices [[i]] <- sample(1:n, size = train.size, replace = FALSE)
#   }
# } else {
#   ## First option: select different indices at each experiment
#   trainIndices = NULL
# }




# dim(counts)
# View(counts)

## Default for quick test without iterating over all cases
permute <- FALSE


if (parameters$compute) {

  train.test.results.all.variables.per.classifier <- list()

  for (classifier in parameters$classifiers) {

    ## List to store all results
    train.test.results.all.variables <- list()
    #    train.test.results.all.PCs <- list()

    message.with.time("Train/test all computations with constant training proportion: ",
                      signif(parameters$trainingProportion, digits = 3) )
    #message.with.time("\tTrain/test, k=", parameters$knn$k, "; classifier=", classifier)



    #### Associate each analysis of real data with a permutation test ####
    for (permute in c(FALSE, TRUE)) {

      # #### Run classifier with all variables (log2-transformed log counts) ####
      # exp.prefix <-
      #   paste(sep = "_", classifier, parameters$recountID , parameters$data.type["log2norm"] , "allvars")
      # if (permute) {
      #   exp.prefix <- paste(sep = "_", exp.prefix, perm.prefix)
      # }# end if permuted class
      #
      # i <- 1
      if(!is.null(self$trainIndices)){
        for( i in 1:length(self$trainIndices)){
          trainIndex <-self$trainIndices[[i]]
          testIndex <- self$testIndices[[i]]


          train.test.results.all.variables[[exp.prefix]] <-
            one.experiment (
              self,
              trainIndex ,
              testIndex ,

              classifier = classifier,
              #variable.type = variable.type,

              file.prefix = exp.prefix,
              permute = permute,
              k = parameters$knn$k,
              verbose = parameters$verbose
            )

        }

      } else {

        stop("/tyou don't have train/test sets to play with classifier ")
      }



      #### Run classifier with all the principal components ####
      #first.pcs <- data.frame(counts)
      first.pcs <- get("log2norm.prcomp.centred.scaled")
      ## define experiment prefix
      exp.prefix <-
        paste(sep = "_", classifier, parameters$recountID , parameters$data.types["prcomp"])
      if (permute) {
        exp.prefix <- paste(sep = "_", exp.prefix, perm.prefix)
      }# end if permuted class

      train.test.results.all.variables[[exp.prefix]] <-
        one.experiment (
          countTable = first.pcs$x,
          classes = classes,
          trainIndices = trainIndices,
          # trainIndex = sample( log2norm.prcomp.centred.scaled$trainIndex),
          # testIndex = sample(log2norm.prcomp.centred.scaled$testIndex),
          data.type = parameters$data.types["prcomp"],
          classifier = classifier,
          variable.type = "all_PCs",
          trainingProportion = parameters$trainingProportion,

          file.prefix = exp.prefix,
          permute = permute,
          k = parameters$knn$k,
          verbose = parameters$verbose
        )


      #### Run classifier with raw counts (no normalization) ####
      ## we looking here to notice the ipmact of normalization into classifiers
      rawCounts1 <- na.omit(as.data.frame(rawCounts$Counts))
      # dim(rawCounts1)

      ## define experiment prefix
      exp.prefix <-
        paste(sep = "_", classifier, parameters$recountID , parameters$data.types["raw"])
      if (permute) {
        exp.prefix <- paste(sep = "_", exp.prefix, perm.prefix)
      }# end if permuted class

      train.test.results.all.variables[[exp.prefix]] <-
        one.experiment (
          countTable = rawCounts1,
          classes = classes,
          trainIndices = trainIndices,
          # trainIndex= sample( rawCounts$trainIndex),
          # testIndex = sample(rawCounts$testIndex),
          data.type = parameters$data.types["raw"],
          classifier = classifier,
          variable.type = "raw",
          trainingProportion = parameters$trainingProportion,

          file.prefix = exp.prefix,
          permute = permute,
          k = parameters$knn$k,
          verbose = parameters$verbose
        )

      #### Run analysis with raw counts, ordered by DEG (edgeR tool), selecting the 75% genes with the smallest DEG p-value ####
      ## we looking here to notice the ipmact of ordered variables by DEG edgeR tool into classifiers
      DEG.object <- get("DEG.edgeR")
      DEG.Counts <- na.omit(as.data.frame(DEG.object$orderedCountTable))
      sig.variables <- round(ncol(DEG.Counts) * 0.75)
      DEG.Counts <- DEG.Counts[, 1:sig.variables]
      ## define experiment prefix
      exp.prefix <-
        paste(sep = "_", classifier, parameters$recountID , parameters$data.types["DEG"])
      if (permute) {
        exp.prefix <- paste(sep = "_", exp.prefix, perm.prefix)
      }# end if permuted class

      train.test.results.all.variables[[exp.prefix]] <-
        one.experiment (
          countTable = DEG.Counts,
          classes = classes,
          trainIndices = trainIndices,
          # trainIndex = sample( DEG.edgeR$trainIndex ) ,
          # testIndex = sample( DEG.edgeR$testIndex),
          data.type = parameters$data.types["DEG"],
          classifier = classifier,
          variable.type = "DEG",
          trainingProportion = parameters$trainingProportion,
          file.prefix = exp.prefix,
          permute = permute,
          k = parameters$knn$k,
          verbose = parameters$verbose
        )



      #### Run classifier with all the variables importance computed by random forest, and then ordered the variables rely on the most 3/4 imporatnce from all varaibles   ####

      v.importance <- get("ordered.countTable.by.importance")
      ## define experiment prefix
      exp.prefix <-
        paste(sep = "_", classifier, parameters$recountID , parameters$data.types["V.importance"], "allvars")
      if (permute) {
        exp.prefix <- paste(sep = "_", exp.prefix, perm.prefix)
      }# end if permuted class

      train.test.results.all.variables[[exp.prefix]] <-
        one.experiment (
          countTable = as.data.frame(ordered.countTable.by.importance),
          classes = classes,
          trainIndices = trainIndices,
          # trainIndex = sample( log2norm.prcomp.centred.scaled$trainIndex),
          # testIndex = sample(log2norm.prcomp.centred.scaled$testIndex),
          data.type = parameters$data.types["V.importance"],
          classifier = classifier,
          variable.type = "all_v.importance",
          trainingProportion = parameters$trainingProportion,

          file.prefix = exp.prefix,
          permute = permute,
          k = parameters$knn$k,
          verbose = parameters$verbose
        )


    } # end loop over permutation

    #### Plotting the Miscalssification Error rate using all diverse data type all variables with KNN classifier? ####
    ErrorRateBoxPlot(experimentList = train.test.results.all.variables,
                     classifier = classifier,
                     data.type = "diverse_data_type",
                     main = paste(sep="",
                                  classifier, ": all variables vs all PCs,", "\n",
                                  parameters$recountID, ", ",
                                  parameters$iterations, " iterations, ","\n",
                                  data.type = "diverse data type"))

    train.test.results.all.variables.per.classifier[[classifier]] <- train.test.results.all.variables

  } # end of loop over classifiers

  # #### Save an image of the results to enable reloading them withouht recomputing everything ####
  # if (parameters$save.image) {
  #   save.image(file = image.file)
  # }

  ##### if compution not required, you can load the image file without any computations ####
  # } else {
  #   # reload previous results if exist
  #   if (file.exists(image.file)) {
  #     message ("Reloading memory image ", image.file)
  #     load(image.file)
  #   } else {
  #     stop("Cannot reload memory image file ", image.file)
  #   }

} # end else if compute statment


###############################################################################################
#### What is better to using all PCs versus all variables with KNN classifier? ####
# ErrorRateBoxPlot(experimentList = train.test.results.all.variables,
#                  classifier = classifier,
#                  main = paste(sep="",
#                               classifier, ": all variables vs all PCs,", "\n",
#                               parameters$recountID, ", ",
#                               parameters$iterations, " iterations, ","\n",
#                               data.type = "diverse data type"))
#
