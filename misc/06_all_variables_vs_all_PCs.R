
##### Impact of normalisation on classifier perfomances ####

## Here we use all variables. we will further test the impact of
## feature selection in separate scripts.

##

## If requested, reset the parameters for all the study cases
## This is used to re-run the analyses on each study case after
## having changed some parameters in the yaml-specific configuration file
if (project.parameters$global$reload.parameters) {
  project.parameters <- yaml.load_file(configFile)
  project.parameters <- initParallelComputing(project.parameters)
  if (exists("studyCases")) {
    for (recountID in names(studyCases)) {
      parameters <- initRecountID(recountID, project.parameters)
      studyCases[[recountID]]$parameters <- parameters
      for (dataSetName in names(studyCases[[recountID]]$datasetsForTest)) {
        studyCases[[recountID]]$datasetsForTest[[dataSetName]]$parameters <- parameters
      }
      rm(parameters)
    }
  }
}



train.test.results.all.variables.per.classifier <- list()

## Run the whole computation if required
## (this can take several hours depending on the number of datasets and classifier methods)

## Loop over recountIDs
## Loop over classifiers
# classifier <- "svm" ## For quick test
# parameters$classifiers[1] ## For quick test
classifier <- "svm"
for (classifier in project.parameters$global$classifiers) {

  # recountID <- "SRP042620"
  for (recountID in selectedRecountIDs) {

    message.with.time("Running train/test with all variables for recountID\t", recountID)

    train.test.results.all.variables.per.classifier[[recountID]] <- list() ## Instantiate an entry per recountID

    ## Get the recountID-specific parameters from the loaded object
    parameters <- studyCases[[recountID]]$parameters



    ## Instantiate variables
    train.test.results.all.variables <- list() ## List to store all results
    dataTypeLabels <- vector() ## Initiate the list of short labels for the plots
    studyCase <- studyCases[[recountID]]

    #### Associate each analysis of real data with a permutation test ####
    permute <- FALSE ## Default for quick test without iterating over all cases
    for (permute in project.parameters$global$permute) {

      #### Data types to analyse ####
      ## If not specified in config file, take all datasetsForTest
      if (is.null(parameters$data.types.to.test)) {
        parameters$data.types.to.test <- names(studyCase$datasetsForTest)
        project.parameters$global$data.types.to.test <- names(studyCase$datasetsForTest)

        ## Temporary (2018-11-01): discard edgeR and DESeq2-sorted datasets.
        ## Actually these are not data types, variable ordering should be treated as a separate variable, not as a separate dataset.
        parameters$data.types.to.test <- grep(pattern = "_sorted", x = parameters$data.types.to.test, invert = TRUE, value = TRUE)
      }

      ## Loop over data types
      data.type <- "TMM_log2_PC" ## For quick test
      for (data.type in parameters$data.types.to.test) {
        message.with.time("\tRunning train/test with all variables",
                          "\n\trecountID: ", recountID,
                          "\n\tData type: ", data.type,
                          "\n\tClassifier: ", classifier,
                          "\n\tpermuted class labels: ", permute)
        dataset <- studyCase$datasetsForTest[[data.type]]
        # class(dataset)
        # summary(dataset)
        # dim(dataset$dataTable)

        dataTypeLabel <- dataset$dataType
        if (permute) {
          dataTypeLabel <- paste(
            dataTypeLabel,
            project.parameters$global$perm.prefix)
        }
        dataTypeLabels <- append(dataTypeLabels, dataTypeLabel)

        # Check if the dataset belongs to the class DataTableWithTrainTestSets
        if (!is(object = dataset, class2 = "DataTableWithTrainTestSets")) {
          stop("The dataset object does not belong to class DataTableWithTrainTestSets. ")
        }

        ## Check if the train and test indices  were properly defined
        if (is.null(dataset$trainTestProperties$trainIndices)
            || is.null(dataset$trainTestProperties$testIndices)) {
          stop("Train/test sets were not properly defined (required to run classifiers). ")
        }

        #### Run classifier with all variables (log2-transformed log counts) ####
        outParam <- outputParameters(dataset, classifier, permute, createDir = TRUE)

        train.test.results.all.variables[[outParam$filePrefix]] <-
          IterateTrainingTesting(
            dataset,
            classifier = classifier,
            permute = permute
          )

      } # End iterations over data types
    } # End iterations over permutation (TRUE / FALSE)


    #### Plotting the Misclassification Error rate using all diverse data type all variables with the selected classifer ####
    outParam <- outputParameters(dataset, classifier = classifier, permute = FALSE, createDir = TRUE)
    outParam$filePrefix <- paste(sep = "_", recountID, classifier, "normalisation_impact")

    ErrorRateBoxPlot(experimentList = train.test.results.all.variables,
                     classifier = classifier,
                     experimentLabels = dataTypeLabels,
                     horizontal = TRUE,
                     fig.height = 6,
                     expMisclassificationRate = dataset$randExpectedMisclassificationRate,
                     # boxplotFile = NULL,
                     boxplotFile = file.path(
                       outParam$resultDir, "figures",
                       paste(sep = "", outParam$filePrefix, ".pdf")),
                     main = paste0(parameters$short_label,
                                   "\n(", parameters$recountID, ")",
                                   "\n", toupper(classifier),
                                   "; ", project.parameters$global$iterations, " iterations")
    )

    ## We store the training-testing result in a single list for further processing
    train.test.results.all.variables.per.classifier[[recountID]][[classifier]] <- train.test.results.all.variables

  } # end loop over RecountIDs
} # end loop over classifiers


## Save the results in a separate object, that can be reloaded later
## Define the path to the memory image for this test (compare classifier whenn they use all variables as features)
save.result.file <- file.path(
  project.parameters$global$dir$memoryImages,
  paste0(
    paste(collapse = "-", selectedRecountIDs),
    "_", project.parameters$global$feature,
    "_normalisation_impact",
    "_", paste(collapse = "-", project.parameters$global$classifiers),
    "_results.Rdata"))
message.with.time(
  "Saving results after eval of normalisation impact: ",
  save.result.file)
dir.create(project.parameters$global$dir$memoryImages, showWarnings = FALSE, recursive = TRUE)
save(train.test.results.all.variables.per.classifier, file = save.result.file)



#### Summarize results for all the study cases and all the classifiers
recountIDs <- names(train.test.results.all.variables.per.classifier)
recountID <- recountIDs[1]
classifiers <- names(train.test.results.all.variables.per.classifier[[recountID]])
# par(mfrow = c(length(recountIDs), length(classifiers)))
for (recountID in recountIDs) {
  for (classifier in classifiers) {
    studyCase <- studyCases[[recountID]]

    message("Summarizing results. RecountID\t", recountID, "\tclassifier\t", classifier)
    experimentList <- train.test.results.all.variables.per.classifier[[recountID]][[classifier]]

    #### Summary table ####
    TTsummary <- SummarizeTrainTestResults(experimentList = experimentList)
    names(TTsummary)
    # expLabels <- TTsummary$experimentLabels
    # expLabels <- sub(pattern = paste0(recountID, "_"), replacement = "", x = expLabels)
    # expLabels <- sub(pattern = paste0(classifier, "_"), replacement = "", x = expLabels)
    # expLabels <- sub(pattern = paste0(project.parameters$global$svm$kernel, "_"), replacement = "", x = expLabels)

    expLabels <- dataTypeLabels
    expLabels <- sub(pattern = "RLE", replacement = "DESeq", x = expLabels)

    #### Error box plots ####

    ## Catch first dataset to get parameters
    parameters <- studyCase$parameters
    stuyCaseLabel <- parameters$short_label
    dataset <- studyCase$datasetsForTest[[1]]
    expMER <- dataset$randExpectedMisclassificationRate
    outParam <- outputParameters(dataset, classifier = classifier, permute = FALSE, createDir = TRUE)
    outParam$filePrefix <- paste(sep = "_", recountID, classifier, "normalisation_impact")
    train.test.results.all.variables <- train.test.results.all.variables.per.classifier[[recountID]][[classifier]]



    ErrorRateBoxPlot(experimentList = train.test.results.all.variables,
                     classifier = classifier,
                     experimentLabels = expLabels,
                     horizontal = TRUE,
                     fig.height = 6,
                     expMisclassificationRate = dataset$randExpectedMisclassificationRate,
                     # boxplotFile = NULL,
                     boxplotFile = file.path(
                       outParam$resultDir, "figures",
                       paste(sep = "", outParam$filePrefix, ".pdf")),
                     main = paste0(parameters$short_label,
                                   "\n(", parameters$recountID, ")",
                                   "\n", toupper(classifier),
                                  "; ", project.parameters$global$iterations, " iterations")
                     )
  }
}
# par(mfrow = c(1, 1))

# stop("HELLO")

## Save a memory image that can be re-loaded next time to avoid re-computing all the normalisation and so on.
if (project.parameters$global$save.image) {
  ## Define the path to the memory image for this test (compare classifier whenn they use all variables as features)
  mem.image.file <- file.path(
    project.parameters$global$dir$memoryImages,
    paste0(
      paste(collapse = "-", selectedRecountIDs),
      "_", project.parameters$global$feature,
      "_normalisation_impact",
      "_", paste(collapse = "-", project.parameters$global$classifiers),
      "_memory-image.Rdata")
  )

  dir.create(project.parameters$global$dir$memoryImages, showWarnings = FALSE, recursive = TRUE)
  message.with.time("Saving memory image after eval of all variables: ", mem.image.file)
  save.image(file = mem.image.file)
}

###############################################################################################
# ErrorRateBoxPlot(experimentList = train.test.results.all.variables,
#                  classifier = classifier,
#                  main = paste(sep="",
#                               classifier, ": all variables vs all PCs,", "\n",
#                               parameters$recountID, ", ",
#                               parameters$iterations, " iterations, ","\n",
#                               data.type = "diverse data type"))
#

message.with.time("Finished script 06_all_variables_vs_all_PCs.R")
