
##### All variables versus all PCs. #####
##
## QUESTION: is it better to use the PCAs-transformed data, and, if so, is it better to use a subset of the first components or all the components ?
## For the time being we test this with only one classifier (KNN, default k)  but we will come back to it with other classifiers later.
## IMPORTANT NOTE : i would like to pay your intention for here we should take " data.type, so that we will not give the user to
## choose the data.type in return we will pass the data.type for each experiment.
## Choice of the classifier

## Define the path to the memory image for this test (compare classifier whenn they use all variables as features)
memory.image.file <- file.path(parameters$dir$memoryImages, "classifier_eval_with_all_variables.Rdata")

## Run the whole computation if required
## (this can take several hours depending on the number of datasets and classifier methods)
if (parameters$compute) {

  train.test.results.all.variables.per.classifier <- list()

  ## Loop over recountIDs
  for (recountID in selectedRecountIDs) {

    train.test.results.all.variables.per.classifier[[recountID]] <- list() ## Instantiate an entry per recountID

    ## Get the recountID-specific parameters from the loaded object
    parameters <- studyCases[[recountID]]$parameters

    #### TEMPORARY FOR DEBUG ####
    # parameters$classifiers <- "svm"
    # parameters$data.types.to.test <- "log2norm"

    message.with.time("Running train/test with all variables for recountID\t", recountID)
    ## Loop over classifiers
    classifier <- "svm" ## For quick test

    for (classifier in parameters$classifiers) {

      ## List to store all results
      train.test.results.all.variables <- list()
      #    train.test.results.all.PCs <- list()

      #message.with.time("\tTrain/test, k=", parameters$knn$k, "; classifier=", classifier)


      #### Associate each analysis of real data with a permutation test ####
      permute <- FALSE ## Default for quick test without iterating over all cases
      for (permute in parameters$permute) {

        ## Loop over data types
        data.type <- "log2norm" ## For test
        for (data.type in parameters$data.types.to.test) {
          message.with.time("\tRunning train/test with all variables",
                            "\n\trecountID: ", recountID,
                            "\n\tClassifier: ", classifier,
                            "\n\tpermuted class labels: ", permute,
                            "\n\tData type: ", data.type)
          dataset <- studyCases[[recountID]][[data.type]]
          # class(dataset)
          # summary(dataset)

          # Check if the dataset belongs to the class countTableWithTrainTestSets
          if (!is(object = dataset, class2 = "countTableWithTrainTestSets")) {

            ## Check if the train and test indices  were properly defined
            if (is.null(dataset$trainTestProperties$trainIndices) || is.null(dataset$trainTestProperties$testIndices)) {
              stop("you don't have train/test sets to play with classifier ")
            }
          }

          #### Run classifier with all variables (log2-transformed log counts) ####
          exp.prefix <-
            paste(sep = "_", recountID, classifier, dataset$dataType)
          if (permute) {
            exp.prefix <- paste(sep = "_", exp.prefix, parameters$perm.prefix)
          }# end if permuted class
          # print(exp.prefix)

          train.test.results.all.variables[[exp.prefix]] <-
            one.experiment (
              dataset,
              classifier = classifier,
              permute = permute#,
              # k = parameters$knn$k,
              # verbose = parameters$verbose
            )

        } # End iterations over dataset
      } # End iterations over permutation


      # #### Run classifier with all the principal components ####
      # #first.pcs <- data.frame(counts)
      # #first.pcs <- get("log2norm.prcomp.centred.scaled")
      # first.pcs <- PCsWithTrainTestSets(studyCases$log2norm)
      #
      # ## define experiment prefix
      # exp.prefix <-
      #   paste(sep = "_", classifier, first.pcs$ID , first.pcs$dataType)
      # if (permute) {
      #   exp.prefix <- paste(sep = "_", exp.prefix, perm.prefix)
      # }# end if permuted class
      #
      # train.test.results.all.variables[[exp.prefix]] <-
      #   one.experiment (
      #     first.pcs,
      #     # classes = classes,
      #     # trainIndices = trainIndices,
      #     # trainIndex = sample( log2norm.prcomp.centred.scaled$trainIndex),
      #     # testIndex = sample(log2norm.prcomp.centred.scaled$testIndex),
      #     # data.type = parameters$data.types["prcomp"],
      #     classifier = classifier,
      #     # variable.type = "all_PCs",
      #     # trainingProportion = parameters$trainingProportion,
      #
      #     file.prefix = exp.prefix,
      #     permute = permute,
      #     k = parameters$knn$k,
      #     verbose = parameters$verbose
      #   )


      #### Run classifier with raw counts (no normalization) ####
      ## we looking here to notice the ipmact of normalization into classifiers
      # dataset2 <- studyCases$filtered
      # # dim(rawCounts1)
      #
      # ## define experiment prefix
      # exp.prefix <-
      #   paste(sep = "_", classifier, parameters$recountID , parameters$data.types["raw"])
      # if (permute) {
      #   exp.prefix <- paste(sep = "_", exp.prefix, perm.prefix)
      # }# end if permuted class
      #
      # train.test.results.all.variables[[exp.prefix]] <-
      #   one.experiment (
      #     dataset2,
      #     # classes = classes,
      #     # trainIndices = trainIndices,
      #     # trainIndex= sample( rawCounts$trainIndex),
      #     # testIndex = sample(rawCounts$testIndex),
      #     # data.type = parameters$data.types["raw"],
      #     classifier = classifier,
      #     # variable.type = "raw",
      #     # trainingProportion = parameters$trainingProportion,
      #
      #     file.prefix = exp.prefix,
      #     permute = permute,
      #     k = parameters$knn$k,
      #     verbose = parameters$verbose
      #   )

      #### Run analysis with raw counts, ordered by DEG (edgeR tool), selecting the 75% genes with the smallest DEG p-value ####
      ## we looking here to notice the ipmact of ordered variables by DEG edgeR tool into classifiers
      # DEG.object <- get("DEG.edgeR")
      # DEG.Counts <- na.omit(as.data.frame(DEG.object$orderedCountTable))
      # sig.variables <- round(ncol(DEG.Counts) * 0.75)
      # DEG.Counts <- DEG.Counts[, 1:sig.variables]
      # ## define experiment prefix
      # exp.prefix <-
      #   paste(sep = "_", classifier, parameters$recountID , parameters$data.types["DEG"])
      # if (permute) {
      #   exp.prefix <- paste(sep = "_", exp.prefix, perm.prefix)
      # }# end if permuted class
      #
      # train.test.results.all.variables[[exp.prefix]] <-
      #   one.experiment (
      #     countTable = DEG.Counts,
      #     classes = classes,
      #     trainIndices = trainIndices,
      #     # trainIndex = sample( DEG.edgeR$trainIndex ) ,
      #     # testIndex = sample( DEG.edgeR$testIndex),
      #     data.type = parameters$data.types["DEG"],
      #     classifier = classifier,
      #     variable.type = "DEG",
      #     trainingProportion = parameters$trainingProportion,
      #     file.prefix = exp.prefix,
      #     permute = permute,
      #     k = parameters$knn$k,
      #     verbose = parameters$verbose
      #   )
      #
      #

      #### Run classifier with all the variables importance computed by random forest, and then ordered the variables rely on the most 3/4 imporatnce from all varaibles   ####

      # v.importance <- get("ordered.countTable.by.importance")
      # ## define experiment prefix
      # exp.prefix <-
      #   paste(sep = "_", classifier, parameters$recountID , parameters$data.types["V.importance"], "allvars")
      # if (permute) {
      #   exp.prefix <- paste(sep = "_", exp.prefix, perm.prefix)
      # }# end if permuted class
      #
      # train.test.results.all.variables[[exp.prefix]] <-
      #   one.experiment (
      #     countTable = as.data.frame(ordered.countTable.by.importance),
      #     classes = classes,
      #     trainIndices = trainIndices,
      #     # trainIndex = sample( log2norm.prcomp.centred.scaled$trainIndex),
      #     # testIndex = sample(log2norm.prcomp.centred.scaled$testIndex),
      #     data.type = parameters$data.types["V.importance"],
      #     classifier = classifier,
      #     variable.type = "all_v.importance",
      #     trainingProportion = parameters$trainingProportion,
      #
      #     file.prefix = exp.prefix,
      #     permute = permute,
      #     k = parameters$knn$k,
      #     verbose = parameters$verbose
      #   )
      #


      #### Plotting the Miscalssification Error rate using all diverse data type all variables with KNN classifier? ####
      ErrorRateBoxPlot(experimentList = train.test.results.all.variables,
                       classifier = classifier,
                       data.type = "diverse-data-types",
                       main = paste(sep="",
                                    parameters$recountID,
                                    "; ", classifier,
                                    "\nall variables; ",
                                    parameters$iterations, " iterations")
                       )
      train.test.results.all.variables.per.classifier[[recountID]][[classifier]] <- train.test.results.all.variables

    } # end loop over classifiers
  } # end loop over recountIDs

  ## Save a memory image that can be re-loaded next time to avoid re-computing all the normalisation and so on.
  if (parameters$save.image) {
    dir.create(parameters$dir$memoryImages, showWarnings = FALSE, recursive = TRUE)
    message.with.time("Saving memory image after eval of all variables: ", memory.image.file)
    save.image(file = memory.image.file)
  }

  ##### if computation not required, you can load the image file without any computations ####
  # } else {
  #   # reload previous results if exist
  #   if (file.exists(memory.image.file)) {
  #     message ("Reloading memory image ", memory.image.file)
  #     load(memory.image.file)
  #   } else {
  #     stop("Cannot reload memory image file ", memory.image.file)
  #   }

}  # end of "if compute"


  #} # end else if compute statment


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
