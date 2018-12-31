## Run training-testing iterations with increasing number of principal components.

##### define the file to store memory Image for " the Number of PCs" test #####
image.dir <- project.parameters$global$dir$memoryImages
dir.create(image.dir, showWarnings = FALSE, recursive = TRUE)

## Define the path to the memory image for this test (compare classifier whenn they use first PCs as features)
featureType <- project.parameters$global$feature
image.file <- file.path(
  project.parameters$global$dir$memoryImages,
  paste(sep = "", "classif_eval_nb-of-PCs_",
        paste(collapse = "-", selectedRecountIDs),
        "_", featureType,
        "_", Sys.Date(), ".Rdata"))

#### Reset the parameters for all the study cases. ####
## This is used to re-run the analyses on each study case after having changed some parameters in the yaml-specific configuration file.
if (project.parameters$global$reload.parameters) {
  project.parameters <- yaml.load_file(configFile)
  project.parameters <- initParallelComputing(project.parameters)
  if (exists("studyCases")) {
    for (recountID in names(studyCases)) {
      parameters <- initRecountID(recountID, project.parameters)
      studyCases[[recountID]]$parameters <- parameters
      #      for (dataSetName in names(studyCases[[recountID]]$datasetsForTest)) {
      # pc.numbers <- c(
      #   2, 3, 4, 5, 6, 7,
      #   seq(from = 10,
      #       to = nrow(studyCases[[recountID]]$datasetsForTest$log2normPCs$dataTable) - 1, by = 10), nrow(studyCases[[recountID]]$datasetsForTest$log2normPCs$dataTable))
      # studyCases[[recountID]]$parameters$pc.numbers <- pc.numbers
      # studyCases[[recountID]]$datasetsForTest[[dataSetName]]$parameters <- parameters
      #      }
      #  print (studyCases[[recountID]]$parameters$dir$tablesDetail)
    }
  }
}

train.test.results.PCs <- list()


#### Iterate over recountIDs ####
classifiers <- project.parameters$global$classifiers
classifier <- "svm" ## For testing
classifier <- "knn" ## For testing
for (classifier in classifiers) {
  ## Instantiate a list to store all results for the current classifier
  train.test.results.PCs[[classifier]] <- list()

  recountID <- selectedRecountIDs[[1]]
  for (recountID in selectedRecountIDs) {
    message("Impact of PCs\t", recountID)

    ## Instantiate a list to store all results for the current recountID
    train.test.results.PCs[[classifier]][[recountID]] <- list()


    #  train.test.results.all.PCs.per.classifier[[recountID]] <- list() ## Instantiate an entry per recountID
    parameters <- studyCases[[recountID]]$parameters
    #studyCases$PCsVar <- PCsWithTrainTestSets(studyCases$filtered)

    #### Iterate over PCs.to.test ####
    data.type <- "TMM_log2_PC"
    for (data.type in project.parameters$global$PCs.to.test) {
      # datasetFile <-  file.path(
      #   project.parameters$global$dir$memoryImages,
      #   paste0("loaded_dataset_", recountID, "_", data.type, ".Rdata"))
      # load(datasetFile, verbose = TRUE)
      # dataset

      #      dataset <- studyCases[[recountID]]$datasetsForTest$log2normPCs$dataTable
      dataset <- studyCases[[recountID]]$datasetsForTest[[data.type]]
      dataTable <- dataset$dataTable
      # dim(dataTable)
      # View(dataTable)

      #### Save original data table to restore it after the tests ####
      original.dataTable <- dataset$dataTable

      #### Define PC numbers depending on dataset size ####
      pc.numbers <- c(
        2, 3, 4, 5, 6, 7,
        seq(from = 10,
            to = nrow(dataTable) - 1, by = 10), nrow(dataTable))
      # length(pc.numbers)

      ## TEMPORARY FOR DEBUGGING
      # pc.numbers <- pc.numbers[1:8]


      studyCases[[recountID]]$parameters$pc.numbers <- pc.numbers
      dataset$parameters$pc.numbers <- pc.numbers
      # studyCases[[recountID]]$datasetsForTest[[dataSetName]]$parameters <- parameters


      #### Iterate over permute FALSE/TRUE ####
      permute <- FALSE
      for (permute in project.parameters$global$permute) {

        #### Iterate over PC numbers ####
        pc.nb <- dataset$parameters$pc.numbers[[1]]
        for (pc.nb in studyCases[[recountID]]$parameters$pc.numbers) {

          message("\n\trecountID:\t", recountID,
                  "\n\tclassifier:\t", classifier,
                  "\n\tdata.type:\t", data.type,
                  "\n\tpermute:\t", permute,
                  "\n\tpc.nb:\t", pc.nb
          )


          ## Select the first N principal components
          first.pcs <- data.frame(t(dataset$prcomp$x[,1:pc.nb]))
          rownames(first.pcs) <- rownames(t(dataset$prcomp$x[,1:pc.nb]))
          dataset$dataTable <- first.pcs

          ## Define experiment prefix
          currentOutputParameters <- outputParameters(dataset = dataset, classifier = classifier, permute = permute, createDir = TRUE)
          exp.prefix <-
            paste(sep = "_", classifier, dataset$ID  , dataset$dataType, "nb_of_PCs", pc.nb)
          if (permute) {
            exp.prefix <- paste(sep = "_", exp.prefix,parameters$perm.prefix)
          }# end if permuted class

          message.with.time("\t", "Experiment prefix: ", exp.prefix)

          ## Define output parameters
          outParam <- outputParameters(dataset, classifier, permute, createDir = TRUE)
          outParam$filePrefix <- paste0(outParam$filePrefix, "_", pc.nb, "_PCs")
          #### Run a training/testing experiment ####
          train.test.results.PCs[[classifier]][[recountID]][[exp.prefix]] <-
            IterateTrainingTesting(
              dataset,
              classifier = classifier,
              file.prefix = outParam$filePrefix,
              permute = permute
            )
        } # end iteration over nb of PCs

      } # end of permutation


      #### Print the results of the effect of the number of PCs on the efficiancy of each classifier classifier ####
      outParam <- outputParameters(
        dataset,
        classifier = paste0(classifier, "_first_PCs"),
        permute = FALSE, createDir = TRUE)
      dir.create(path = file.path(outParam$resultDir, "figures"), showWarnings = FALSE, recursive = FALSE)
      ErrorRateBoxPlot(experimentList = train.test.results.PCs[[classifier]][[recountID]],
                       classifier = classifier,
                       horizontal = TRUE,
                       main = paste0(
                         "PC selection; ", parameters$recountID, "; ", dataset$dataType, "\n",
                         classifier, "; ", parameters$iterations , "iterations,"),
                       boxplotFile = file.path(
                         outParam$resultDir, "figures",
                         paste(sep = "", outParam$filePrefix, ".pdf"))
      )


      #  train.test.results.all.PCs.per.classifier[[recountID]][[classifier]] <- train.test.results.PCs[[classifier]][[recountID]]
    }  # end of loop over classifiers

    ## Restore original data table
    dataset$dataTable <- original.dataTable
  } # end loop over datasets
} # end loop over recountIDs

#### Summarize results for all recountIDs ####
for (classifier in classifiers) {
  allRecountIDsSummary <- list()
  allRecountIDsSummary$perRecountID <- list()
  allRecountIDsSummary$summary <- data.frame()
  # rownames(allRecountIDsSummary$summary) <- names(train.test.results.PCs)
  recountIDs <- names(train.test.results.PCs)
  for (i in 1:length(recountIDs)) {
    recountID <- recountIDs[i]
    oneRecountIDSummary <- SummarizeTrainTestResults(experimentList = train.test.results.PCs[[classifier]][[recountID]])
    allRecountIDSummary$perRecountID[[recountID]] <- oneRecountIDSummary

    allRecountIDsSummary$summary <- rbind(
      allRecountIDsSummary$summary,
      data.frame(
        test.min = min(oneRecountIDSummary$testing.MER.summary$mean),
        test.which.min = which.min(oneRecountIDSummary$testing.MER.summary$mean),
        test.which.min.name = rownames(oneRecountIDSummary$testing.MER.summary)[which.min(oneRecountIDSummary$testing.MER.summary$mean)],
        test.max = max(oneRecountIDSummary$testing.MER.summary$mean),
        test.which.max = which.max(oneRecountIDSummary$testing.MER.summary$mean),
        test.which.max.name = rownames(oneRecountIDSummary$testing.MER.summary)[which.max(oneRecountIDSummary$testing.MER.summary$mean)]
      )
    )
  }
  rownames(allRecountIDsSummary$summary) <- recountIDs
}

## Save the results in a separate object, that can be reloaded later
## Define the path to the memory image for this test (compare classifier whenn they use all variables as features)
dir.create(project.parameters$global$dir$memoryImages, showWarnings = FALSE, recursive = TRUE)
featureType <- project.parameters$global$feature
save.result.file <- file.path(
  project.parameters$global$dir$memoryImages,
  paste0(
    "feature_selection_first_PCs",
    "_", paste(collapse = "-", selectedRecountIDs),
    "_", featureType,
    ".Rdata"))
save(train.test.results.PCs, file = save.result.file)
message.with.time(
  "Saved results after eval of feature selection with first PCs: ",
  save.result.file)


message.with.time("Finished script 07_PCA_impact_of_PC_number.R")

