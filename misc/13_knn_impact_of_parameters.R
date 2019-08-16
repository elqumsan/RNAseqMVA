#### Test the impact of k on KNN performances ####

## List to store the results for all recountIDs
train.test.results.knn.k.values <- list()

#### Loop over recountIDs ####
recountID <- names(studyCases)[1]
for (recountID in selectedRecountIDs) {
  message.with.time("Running train/test with all variables to test imapct of knn's parameters for recountID\t", recountID)

  ## Initialize a list to store all results for the current recountID
  train.test.results.all.variables.knn <- list()

  ## Get the recountID-specific parameters from the loaded object
  parameters <- studyCases[[recountID]]$parameters

  #### Check the data types to analyse ####
  ## If not specified in config file, take all datasetsForTest
  if (is.null(parameters$data.types.to.test)) {
    parameters$data.types.to.test <- names(studyCase$datasetsForTest)

    ## Temporary (2018-11-01): discard edgeR and DESeq2-sorted datasets.
    ## Actually these are not data types, variable ordering should be treated as a separate variable, not as a separate dataset.
    parameters$data.types.to.test <- grep(pattern = "DESeq2", x = parameters$data.types.to.test, invert = TRUE, value = TRUE)
    parameters$data.types.to.test <- grep(pattern = "edgeR", x = parameters$data.types.to.test, invert = TRUE, value = TRUE)

    ## Update global parameters (required for the comparative plots)
    project.parameters$global$data.types.to.test <- parameters$data.types.to.test

  }


  ## Make sure that the classifier variable is KNN
  classifier <- "knn"

  #### Loop over label permutation ####
  permute <- FALSE ## Default for quick test without iterating over all cases
  for (permute in project.parameters$global$permute) {


    #### Loop over values of k ####
    knn.k <- 3
    for (knn.k in project.parameters$global$knn$k_values) {

      #### Loop over data types ####
      data.type <- "TMM_log2" ## For test
      for (data.type in project.parameters$global$data.types.to.test) {
        message.with.time("\tRunning train/test with all variables",
                          "\n\trecountID:\t", recountID,
                          "\n\tClassifier:\t", classifier,
                          "\n\tpermuted:\t", permute,
                          "\n\tk:        \t", knn.k,
                          "\n\tData type:\t", data.type)
        dataset <- studyCases[[recountID]]$datasetsForTest[[data.type]]
        dataset$parameters$knn$k <- knn.k
        # class(dataset)
        # summary(dataset)

        # Check if the dataset belongs to the class DataTableWithTrainTestSets
        if (!is(object = dataset, class2 = "DataTableWithTrainTestSets")) {

          ## Check if the train and test indices  were properly defined
          if (is.null(dataset$trainTestProperties$trainIndices) || is.null(dataset$trainTestProperties$testIndices)) {
            stop("you don't have train/test sets to play with classifier ")
          }
        }

        ## Define output parameters
        outParam <- outputParameters(dataset, classifier, permute, createDir = TRUE)
        exp.prefix <- outParam$filePrefix
        ## exp.prefix <- paste0(recountID, "_", data.type, "_knn_k", knn.k)

        #### Run a training/testing experiment ####
        train.test.results.all.variables.knn[[exp.prefix]] <-
          IterateTrainingTesting(
            dataset,
            classifier = classifier,
            file.prefix = outParam$filePrefix,
            permute = permute)

      } # End iterations over dataset

    } # End iterations over permutation

  } # end iteration over k values

  #### Plotting the Misclassification Error rate using all diverse data type all variables with KNN classifier? ####
  outParam <- outputParameters(
    dataset,
    classifier = "knn impact of k",
    permute = FALSE, createDir = TRUE)
  dir.create(path = file.path(outParam$resultDir, "figures"), showWarnings = FALSE, recursive = FALSE)
  ErrorRateBoxPlot(experimentList = train.test.results.all.variables.knn,
                   classifier = classifier,
                   horizontal = TRUE,
                   boxplotFile = file.path(
                     outParam$resultDir, "figures",
                     paste(sep = "", outParam$filePrefix, ".pdf")),
                   main = paste(sep = "",
                                parameters$recountID,
                                "; ", classifier, "; ",
                                "\nall variables; ",
                                project.parameters$global$iterations, " iterations")
  )
  train.test.results.knn.k.values[[recountID]] <- train.test.results.all.variables.knn

} # end loop over recountIDs

for (recountID in names(studyCases)) {
  outParam <- outputParameters(
    dataset,
    classifier = "knn impact of k",
    permute = FALSE, createDir = TRUE)
  ErrorRateBoxPlot(experimentList = train.test.results.all.variables.knn[[recountID]],
                   classifier = classifier,
                   horizontal = TRUE,
                   boxplotFile = file.path(
                     outParam$resultDir, "figures",
                     paste(sep = "", outParam$filePrefix, ".pdf")),
                   main = paste(sep = "",
                                parameters$recountID,
                                "; ", classifier, "; ",
                                "\nall variables; ",
                                project.parameters$global$iterations, " iterations"))
}


## Save the results in a separate object, that can be reloaded later
## Define the path to the memory image for this test (compare classifier whenn they use all variables as features)
save.result.file <- file.path(project.parameters$global$dir$memoryImages, "knn_impact_of_k_result.Rdata")
dir.create(project.parameters$global$dir$memoryImages, showWarnings = FALSE, recursive = TRUE)
save(train.test.results.knn.k.values, file = save.result.file)
message.with.time("Saving results  after eval of impact of k on knn: ", save.result.file)

## Save a memory image that can be re-loaded next time to avoid re-computing all the normalisation and so on.
if (project.parameters$global$save.image) {
  message.with.time("Saving memory image after eval of impact of k on knn: ", memory.image.file)
  save.image(file = memory.image.file)
}



message.with.time("Finished script 13_knn_impact_of_parameters.R")
