
##### Test the impact of k on KNN performances

## Define the path to the memory image for this test (compare classifier whenn they use all variables as features)
memory.image.file <- file.path(project.parameters$global$dir$memoryImages, "knn_impact_of_k.Rdata")

## For debug: reset the parameteres for all the study cases
## This is used to re-run the analyses on each study case after having changed some parameters in the yaml-specific configuration file
project.parameters <- yaml.load_file(configFile)
project.parameters <- initParallelComputing(project.parameters)
if (exists("studyCases")) {
  message("\tResetting parameters for the ", length(studyCases), " study cases")
  recountID <- names(studyCases)[1]
  for (recountID in names(studyCases)) {
    parameters <- initRecountID(recountID, project.parameters)
    studyCases[[recountID]]$parameters <- parameters
    for (dataSetName in names(studyCases[[recountID]]$datasetsForTest)) {
      studyCases[[recountID]]$datasetsForTest[[dataSetName]]$parameters <- parameters
    }
    #  print (studyCases[[recountID]]$parameters$dir$tablesDetail)
  }
}

## List to store the results for all recountIDs
train.test.results.knn.k.values <- list()

#### Loop over recountIDs ####
recountID <- names(studyCases)[1]
for (recountID in selectedRecountIDs) {

  ## List to store all results for the current recountID
  train.test.results.all.variables.knn <- list()

  ## Get the recountID-specific parameters from the loaded object
  parameters <- studyCases[[recountID]]$parameters

  message.with.time("Running train/test with all variables to test imapct of knn's parameters for recountID\t", recountID)
  ## Loop over classifiers
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
