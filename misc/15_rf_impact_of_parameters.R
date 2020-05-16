##### Impact of ntree, mtry and nodesize parameters for the randomForest classifier #####
##
## QUESTION: is it better to amendment these rf's aparameters in order to optimize the rf effeciency, and, if so, what is the best values for those parameter in order to increase the leverage of randomForest classifier. ?
## For the time being we test the rf with default values for ntree = 500, mtry =sqrt(number of features = 55000) and nodesize = 1
## IMPORTANT NOTE : i would like to pay your intention for here we should take " data.type, so that we will not give the user to
## choose the data.type in return we will pass the data.type for each experiment.
## we wil test impact of the ntree = 1000
##                           mtry = 500
##                           nodesize = 20000

## Define the path to the memory image for this test (compare classifier whenn they use all variables as features)
memory.image.file <- file.path(project.parameters$global$dir$memoryImages, "rf_impact_of_parameters.Rdata")

## For debug: reset the parameteres for all the study cases
## This is used to re-run the analyses on each study case after having changed some parameters in the yaml-specific configuration file
project.parameters <- yaml.load_file(configFile)
# project.parameters <- initParallelComputing(project.parameters)
if (exists("studyCases")) {
  #  for (recountID in names(studyCases)) {
  parameters <- initRecountID(recountID, project.parameters)
  studyCases[[recountID]]$parameters <- parameters
  for (dataSetName in names(studyCases[[recountID]]$datasetsForTest)) {
    studyCases[[recountID]]$datasetsForTest[[dataSetName]]$parameters <- parameters
  }
  #  print (studyCases[[recountID]]$parameters$dir$tablesDetail)
  #  }
}



## Run the whole computation if required
## (this can take several hours depending on the number of datasets and classifier methods)
train.test.results.all.variables.per.rf <- list()

train.test.results.all.variables.per.rf[[recountID]] <- list() ## Instantiate an entry per recountID

## Get the recountID-specific parameters from the loaded object
parameters <- studyCases[[recountID]]$parameters

#### TEMPORARY FOR DEBUG ####
# parameters$classifiers <- "rf"
# parameters$data.types.to.test <- "log2norm"
# parameters$data.types.to.test <- "log2normPCs"

message.with.time("Running train/test with all variables to test imapct of rf's parameters for recountID\t", recountID)
## Loop over classifiers
classifier <- "rf"

for (mtry in project.parameters$global$rf$mtry) {

  for (ntree in project.parameters$global$rf$ntree) {

    ## List to store all results
    train.test.results.all.variables.rf <- list()
    #    train.test.results.all.PCs <- list()

    #message.with.time("\tTrain/test, k=", parameters$knn$k, "; classifier=", classifier)


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
        dataset$parameters$rt$ntree <- ntree
        dataset$parameters$rf$mtry <- mtry
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
        outputParam  <-outputParameters(dataset , classifier )
        outputParam$filePrefix <- paste(outputParam$filePrefix, ntree,"ntree",mtry,"mtry" ,sep = "_")
        # exp.prefix <- filePrefix(dataset,classifier, permute)
        #  paste(sep = "_", recountID, classifier, dataset$dataType)
        if (permute) {
          #  exp.prefix <- paste(sep = "_", exp.prefix, project.parameters$global$perm.prefix)
          outputParam <- outputParameters(dataset,classifier, permute)
          outputParam$filePrefix <- paste(outputParam$filePrefix, ntree,"ntree",mtry,"mtry" ,sep = "_")
        }
        # exp.prefix <-
        #   paste(sep = "_", recountID, classifier, dataset$dataType, svm.kernel)
        # if (permute) {
        #   exp.prefix <- paste(sep = "_", exp.prefix, project.parameters$global$perm.prefix)
        # }# end if permuted class
        # # print(exp.prefix)

        train.test.results.all.variables.rf[[outputParam$filePrefix]] <-
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
    ErrorRateBoxPlot(experimentList = train.test.results.all.variables.rf,
                     classifier = classifier,
                     main = paste(sep = "",
                                  parameters$recountID,
                                  "; ", classifier,
                                  "\nall variables; ",
                                  project.parameters$global$iterations, " iterations")
    )
    train.test.results.all.variables.per.rf[[recountID]][[ntree]] <- train.test.results.all.variables.rf

  } # end iteration over RF-ntree-parameter



} # end of iteration over mtry



## Save a memory image that can be re-loaded next time to avoid re-computing all the normalisation and so on.
if (project.parameters$global$save.image) {
  dir.create(project.parameters$global$dir$memoryImages, showWarnings = FALSE, recursive = TRUE)
  message.with.time("Saving memory image after eval of all variables: ", memory.image.file)
  save.image(file = memory.image.file)
}

