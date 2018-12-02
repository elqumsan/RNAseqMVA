
##### All variables versus all PCs. #####
##
## QUESTION: is it better to use the PCAs-transformed data, and, if so, is it better to use a subset of the first components or all the components ?
## For the time being we test this with only one classifier (KNN, default k)  but we will come back to it with other classifiers later.
## IMPORTANT NOTE : i would like to pay your intention for here we should take " data.type, so that we will not give the user to
## choose the data.type in return we will pass the data.type for each experiment.
## Choice of the classifier

## Define the path to the memory image for this test (compare classifier whenn they use all variables as features)
allVariables.mem.image <- file.path(
  project.parameters$global$dir$memoryImages,
  paste(sep = "", "classif_eval_all_variables_",
      paste(collapse = "-", selectedRecountIDs),
      "_", Sys.Date(), ".Rdata"))


# ## TEMP
reload.mem.image <- TRUE
if (reload.mem.image) {
    studyCases.mem.image <- "~/RNAseqMVA_workspace/memory_images/loaded_studyCases_SRP042620-SRP057196-SRP056295-SRP035988-SRP061240-SRP062966-SRP066834_2018-11-01.Rdata"
    message("Reloading study cases from previously stored memory image")
    message("\t", studyCases.mem.image)
    load(studyCases.mem.image)
}

## For debug: reset the parameters for all the study cases
## This is used to re-run the analyses on each study case after
## having changed some parameters in the yaml-specific configuration file
reload.parameters <- TRUE
if (reload.parameters) {
  project.parameters <- yaml.load_file(configFile)
  project.parameters <- initParallelComputing(project.parameters)
  if (exists("studyCases")) {
    for (recountID in names(studyCases)) {
      parameters <- initRecountID(recountID, project.parameters)
      studyCases[[recountID]]$parameters <- parameters
      for (dataSetName in names(studyCases[[recountID]]$datasetsForTest)) {
        studyCases[[recountID]]$datasetsForTest[[dataSetName]]$parameters <- parameters
      }
      #  print (studyCases[[recountID]]$parameters$dir$tablesDetail)
    }
  }
}

train.test.results.all.variables.per.classifier <- list()

## Run the whole computation if required
## (this can take several hours depending on the number of datasets and classifier methods)
if (project.parameters$global$compute) {


  ## Loop over recountIDs
  # recountID <- "SRP042620"
  for (recountID in selectedRecountIDs) {

    message.with.time("Running train/test with all variables for recountID\t", recountID)

    train.test.results.all.variables.per.classifier[[recountID]] <- list() ## Instantiate an entry per recountID

    ## Get the recountID-specific parameters from the loaded object
    parameters <- studyCases[[recountID]]$parameters


    studyCase <- studyCases[[recountID]]
    ## Loop over classifiers
    # classifier <- "svm" ## For quick test
    for (classifier in parameters$classifiers) {
      short.labels <- vector() ## Initiate the list of short labels for the plots

      ## List to store all results
      train.test.results.all.variables <- list()

      #### Associate each analysis of real data with a permutation test ####
      permute <- FALSE ## Default for quick test without iterating over all cases
      for (permute in project.parameters$global$permute) {

        #### Data types to analyse ####
        ## If not specified in config file, take all datasetsForTest
        if (is.null(parameters$data.types.to.test)) {
          parameters$data.types.to.test <- names(studyCase$datasetsForTest)

          ## Temporary (2018-11-01): discard edgeR and DESeq2-sorted datasets.
          ## Actually these are not data types, variable ordering should be treated as a separate variable, not as a separate dataset.
          parameters$data.types.to.test <- grep(pattern = "_sorted", x = parameters$data.types.to.test, invert = TRUE, value = TRUE)
        }

        ## Loop over data types
        data.type <- "q0.75_log2_PC" ## For quick test
        for (data.type in parameters$data.types.to.test) {
          message.with.time("\tRunning train/test with all variables",
                            "\n\trecountID: ", recountID,
                            "\n\tData type: ", data.type,
                            "\n\tClassifier: ", classifier,
                            "\n\tpermuted class labels: ", permute)
          dataset <- studyCase$datasetsForTest[[data.type]]
          # class(dataset)
          # summary(dataset)

          short.label <- dataset$dataType
          if (permute) {
            short.label <- paste(
              short.label,
              project.parameters$global$perm.prefix)
          }
          short.labels <- append(short.labels, short.label)

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
      } # End iterations over permutation



      #### Plotting the Misclassification Error rate using all diverse data type all variables with KNN classifier? ####
      outParam <- outputParameters(dataset, classifier = classifier, permute = FALSE, createDir = TRUE)
      ErrorRateBoxPlot(experimentList = train.test.results.all.variables,
                       classifier = classifier,
                       experimentLabels = short.labels,
                       horizontal = FALSE,
                       fig.height = 8,
                       expMisclassificationRate = dataset$randExpectedMisclassificationRate,
                       # boxplotFile = NULL,
                       boxplotFile = file.path(
                         outParam$resultDir, "figures",
                         paste(sep = "", outParam$filePrefix, ".pdf")),
                       main = paste(sep = "",
                                    parameters$recountID,
                                    "; ", classifier,
                                    "\n all features; ",
                                    project.parameters$global$iterations, " iterations"))

      ## We store the training-testing result in a single list for further processing
      train.test.results.all.variables.per.classifier[[recountID]][[classifier]] <- train.test.results.all.variables

    } # end loop over classifiers
  } # end loop over RecountIDs
} # end if compute

## Save a memory image that can be re-loaded next time to avoid re-computing all the normalisation and so on.
if (project.parameters$global$save.image) {
  dir.create(project.parameters$global$dir$memoryImages, showWarnings = FALSE, recursive = TRUE)
  message.with.time("Saving memory image after eval of all variables: ", allVariables.mem.image)
  save.image(file = allVariables.mem.image)
}

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

message.with.time("Finished script 06_all_vaariables_vs_all_PCs.R")
