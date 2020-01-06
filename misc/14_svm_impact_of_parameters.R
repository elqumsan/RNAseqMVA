#### Test the impact of the kernel on SVM performances ####

## List to store the results for all recountIDs
train.test.results.all.variables.per.svm <- list()

#### Loop over recountIDs ####
recountID <- names(studyCases)[1]
for (recountID in selectedRecountIDs) {
  message.with.time("Assessing kernel impact on SVM performances\t", recountID)

  ## Initialise a list to store all results for the current recountID
  train.test.results.all.variables.svm <- list()

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


  ## Set the parameters
  classifier <- "svm" ## Make sure that the classiifer is SVM for this script

  #### Loop over label permutation options ####
  permute <- FALSE ## Default for quick test without iterating over all cases
  for (permute in parameters$permute) {

    #### Loop over data types ####
    data.type <- "TMM_log2" ## For test
    for (data.type in parameters$data.types.to.test) {

      #### Loop over kernels ####
      svm.kernel <- "linear" ## Default kernel value
      for (svm.kernel in parameters$svm$kernel_values) {

        ## Verbosity
        message.with.time("\tRunning train/test with all variables",
                          "\n\trecountID: ", recountID,
                          "\n\tClassifier: ", classifier,
                          "\n\tpermuted class labels: ", permute,
                          "\n\tData type: ", data.type)
        dataset <- studyCases[[recountID]]$datasetsForTest[[data.type]]
        dataset$parameters$svm$kernel <- svm.kernel
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
        ## exp.prefix <- paste0(recountID, "_", data.type, "_svm_kernel_", svm.kernel)

        #### Run a training/testing experiment ####
        train.test.results.all.variables.svm[[exp.prefix]] <-
          IterateTrainingTesting(
            dataset,
            classifier = classifier,
            file.prefix = outParam$filePrefix,
            permute = permute)

      } # End iterations over dataset

    } # End iterations over permutation

  } # end iteration over svm_kernels

  #### Plotting the Misclassification Error rate using all diverse data type all variables with KNN classifier? ####
  outParam <- outputParameters(dataset, classifier = "svm kernel comparison", permute = FALSE, createDir = TRUE)
  dir.create(path = file.path(outParam$resultDir, "figures"), showWarnings = FALSE, recursive = FALSE)

  experimentLabels <- names(train.test.results.all.variables.svm)
  experimentLabels <- sub(pattern = paste0(parameters$recountID, "_", parameters$feature, "_svm_"),
                          replacement = "",
                          x = experimentLabels)

  # ## TEMPORARY FIX FOR THE ORDER OF THE EXPERIMENTS
  # kernel <- sub(pattern = "_.*", replacement = "", x = experimentLabels, perl = TRUE)
  # dataType <- sub(pattern = "[^_]+_", replacement = "", x = experimentLabels, perl = TRUE)
  # #  nonPermuted <- grep(pattern = "permLabels", x = experimentLabels, invert = TRUE)
  # #  experimentOrder <- order(dataType[nonPermuted])
  # experimentOrder <- order(dataType)

  ## Replace underscore by space in the labels
  experimentLabels <- sub(pattern = "_", replacement = " ", x = experimentLabels)

  ErrorRateBoxPlot(experimentList = train.test.results.all.variables.svm,
                   classifier = classifier,
                   horizontal = TRUE,
                   experimentLabels = experimentLabels,
                   boxplotFile = file.path(
                     outParam$resultDir, "figures",
                     paste0(parameters$recountID,
                            "_", parameters$feature,
                            "_svm_kernel_comparison",
                            ".pdf")),
                   main = paste(sep = "",
                                parameters$short_label, " (", parameters$recountID, ") ",
                                parameters$feature,
                                "\nImpact of SVM kernel; ",
                                project.parameters$global$iterations, " iterations")
  )
  train.test.results.all.variables.per.svm[[recountID]] <- train.test.results.all.variables.svm

} # end loop over recountIDs


## Save the results in a separate object, that can be reloaded later
## Define the path to the memory image for this test (compare classifier whenn they use all variables as features)
save.result.file <- file.path(
  project.parameters$global$dir$memoryImages,
  paste0(
    paste(collapse = "-", selectedRecountIDs),
    "_", featureType,
    "_svm_impact_of_parameters_result",
    "_", Sys.Date(), "_results.Rdata"))
dir.create(project.parameters$global$dir$memoryImages, showWarnings = FALSE, recursive = TRUE)
message.with.time("Saving results  after eval of SVM kernels: ", save.result.file)
save(train.test.results.all.variables.per.svm, file = save.result.file)

## Save a memory image that can be re-loaded next time to avoid re-computing all the normalisation and so on.
if (project.parameters$global$save.image) {
  memory.image.file <- file.path(
    project.parameters$global$dir$memoryImages,
    paste0(
      paste(collapse = "-", selectedRecountIDs),
      "_", featureType,
      "_", "svm_impact_of_parameters_result",
      "_", Sys.Date(), "_memory-image.Rdata"))
  message.with.time("Saving memory image after eval of SVM kernels: ", memory.image.file)
  save.image(file = memory.image.file)
}



message.with.time("Finished script 14_svm_impact_of_parameters.R")
