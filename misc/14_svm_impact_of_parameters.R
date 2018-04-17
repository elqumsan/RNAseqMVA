
##### All variables versus all PCs. #####
##
## QUESTION: is it better to use the PCAs-transformed data, and, if so, is it better to use a subset of the first components or all the components ?
## For the time being we test this with only one classifier (KNN, default k)  but we will come back to it with other classifiers later.
## IMPORTANT NOTE : i would like to pay your intention for here we should take " data.type, so that we will not give the user to
## choose the data.type in return we will pass the data.type for each experiment.
## Choice of the classifier

## Define the path to the memory image for this test (compare classifier whenn they use all variables as features)
memory.image.file <- file.path(project.parameters$global$dir$memoryImages, "svm_impact_of_parameters.Rdata")


## For debug: reset the parameteres for all the study cases
## This is used to re-run the analyses on each study case after having changed some parameters in the yaml-specific configuration file
project.parameters <- yaml.load_file(configFile)
project.parameters <- initParallelComputing(project.parameters)
if(exists("studyCases")) {
  for (recountID in names(studyCases)) {
    parameters <- initRecountID(recountID, project.parameters)
    studyCases[[recountID]]$parameters <- parameters
    for (dataSetName in names(studyCases[[recountID]]$datasetsForTest)) {
      studyCases[[recountID]]$datasetsForTest[[dataSetName]]$parameters <- parameters
    }
    #  print (studyCases[[recountID]]$parameters$dir$tablesDetail)
  }
}


## Run the whole computation if required
## (this can take several hours depending on the number of datasets and classifier methods)
if (project.parameters$global$compute) {

  train.test.results.all.variables.per.svm <- list()

  ## Loop over recountIDs
  for (recountID in selectedRecountIDs) {

    train.test.results.all.variables.per.svm[[recountID]] <- list() ## Instantiate an entry per recountID

    ## Get the recountID-specific parameters from the loaded object
    parameters <- studyCases[[recountID]]$parameters

    #### TEMPORARY FOR DEBUG ####
    # parameters$classifiers <- "svm"
    # parameters$data.types.to.test <- "log2norm"
    # parameters$data.types.to.test <- "log2normPCs"

    message.with.time("Running train/test with all variables to test imapct of svm's parameters for recountID\t", recountID)
    ## Loop over classifiers
    classifier <- "svm"

    for (svm.kernel in project.parameters$global$svm$kernel) {

      ## List to store all results
      train.test.results.all.variables.svm <- list()
      #    train.test.results.all.PCs <- list()

      #message.with.time("\tTrain/test, k=", parameters$knn$k, "; classifier=", classifier)


      #### Associate each analysis of real data with a permutation test ####
      permute <- FALSE ## Default for quick test without iterating over all cases
      for (permute in project.parameters$global$permute) {

        ## Loop over data types
        data.type <- "log2normPCs" ## For test
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

          #### Run classifier with all variables (log2-transformed log counts) ####
          exp.prefix <-
            paste(sep = "_", recountID, classifier, dataset$dataType, svm.kernel)
          if (permute) {
            exp.prefix <- paste(sep = "_", exp.prefix, project.parameters$global$perm.prefix)
          }# end if permuted class
          # print(exp.prefix)

          train.test.results.all.variables.svm[[exp.prefix]] <-
            IterateTrainingTesting (
              dataset,
              classifier = classifier,
              permute = permute#,
              # k = parameters$knn$k,
              # verbose = parameters$verbose
            )

        } # End iterations over dataset
      } # End iterations over permutation




      #### Plotting the Miscalssification Error rate using all diverse data type all variables with KNN classifier? ####
      ErrorRateBoxPlot(experimentList = train.test.results.all.variables.svm,
                       classifier = classifier,
                       data.type = "diverse-data-types",
                       main = paste(sep="",
                                    parameters$recountID,
                                    "; ", classifier,
                                    "\nall variables; ",
                                    project.parameters$global$iterations, " iterations")
      )
      train.test.results.all.variables.per.svm[[recountID]][[classifier]] <- train.test.results.all.variables.svm

    } # end loop over classifiers
  } # end loop over recountIDs

  ## Save a memory image that can be re-loaded next time to avoid re-computing all the normalisation and so on.
  if (project.parameters$global$save.image) {
    dir.create(project.parameters$global$dir$memoryImages, showWarnings = FALSE, recursive = TRUE)
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

message.with.time("Finished script 13_svm_impact_of_parameters.R")
