
##### Test the impact of kernel on SVM performances

## Define the path to the memory image for this test (compare classifier whenn they use all variables as features)
memory.image.file <- file.path(project.parameters$global$dir$memoryImages, "svm_impact_of_parameters.Rdata")

## For debug: reset the parameteres for all the study cases
## This is used to re-run the analyses on each study case after having changed some parameters in the yaml-specific configuration file
project.parameters <- yaml.load_file(configFile)
project.parameters <- initParallelComputing(project.parameters)
if (exists("studyCases")) {
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
train.test.results.all.variables.per.svm <- list()

## Loop over recountIDs
recountID <- names(studyCases)[1]
for (recountID in selectedRecountIDs) {

  ## List to store all results for the current recountID
  train.test.results.all.variables.svm <- list()

  ## Get the recountID-specific parameters from the loaded object
  parameters <- studyCases[[recountID]]$parameters

  message.with.time("Running train/test with all variables to test imapct of svm's parameters for recountID\t", recountID)
  ## Loop over classifiers
  classifier <- "svm"
  svm.kernel <- "linear"
  for (svm.kernel in project.parameters$global$svm$kernel_values) {


    #### Associate each analysis of real data with a permutation test ####
    permute <- FALSE ## Default for quick test without iterating over all cases
    for (permute in project.parameters$global$permute) {

      ## Loop over data types
      data.type <- "TMM_log2" ## For test
      for (data.type in project.parameters$global$data.types.to.test) {
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


        #### Run classifier with all variables (log2-transformed log counts) ####
        exp.prefix <- outParam$filePrefix

        train.test.results.all.variables.svm[[exp.prefix]] <-
          IterateTrainingTesting(
            dataset,
            classifier = classifier,
            permute = permute)

      } # End iterations over dataset

    } # End iterations over permutation

  } # end iteration over svm_kernels

  #### Plotting the Miscalssification Error rate using all diverse data type all variables with KNN classifier? ####
  ErrorRateBoxPlot(experimentList = train.test.results.all.variables.svm,
                   classifier = classifier,
                   data.type = "diverse-data-types",
                   main = paste(sep = "",
                                parameters$recountID,
                                "; ", classifier, "; ", svm.kernel, " kernel; ",
                                "\nall variables; ",
                                project.parameters$global$iterations, " iterations")
  )
  train.test.results.all.variables.per.svm[[recountID]][[svm.kernel]] <- train.test.results.all.variables.svm

} # end loop over recountIDs


## Save a memory image that can be re-loaded next time to avoid re-computing all the normalisation and so on.
if (project.parameters$global$save.image) {
  dir.create(project.parameters$global$dir$memoryImages, showWarnings = FALSE, recursive = TRUE)
  message.with.time("Saving memory image after eval of all variables: ", memory.image.file)
  save.image(file = memory.image.file)
}



message.with.time("Finished script 13_svm_impact_of_parameters.R")
