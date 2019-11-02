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



train.test.results.PCs <- list()


#### Iterate over recountIDs ####
classifiers <- project.parameters$global$classifiers
# classifier <- "svm" ## For testing
# classifier <- "knn" ## For testing
classifier <- classifiers[1] ## for testing
for (classifier in classifiers) {
  ## Instantiate a list to store all results for the current classifier
  train.test.results.PCs[[classifier]] <- list()

  recountID <- selectedRecountIDs[[1]]
  for (recountID in selectedRecountIDs) {
    #  train.test.results.all.PCs.per.classifier[[recountID]] <- list() ## Instantiate an entry per recountID
    parameters <- studyCases[[recountID]]$parameters

    ## Instantiate a list to store all results for the current recountID
    train.test.results.PCs[[classifier]][[recountID]] <- list()

    #### Iterate over PCs.to.test ####
    data.type <- "TMM_log2_PC"
    for (data.type in project.parameters$global$PCs.to.test) {
      message.with.time("Impact of PCs\t", recountID, " ", parameters$feature, "\t", data.type)
      dataset <- studyCases[[recountID]]$datasetsForTest[[data.type]]

      # datasetFile <-  file.path(
      #   project.parameters$global$dir$memoryImages,
      #   paste0("loaded_dataset_", recountID, "_", data.type, ".Rdata"))
      # load(datasetFile, verbose = TRUE)
      # dataset


      #### Save original data table to restore it after the tests ####
      original.dataTable <- dataset$dataTable
      # dim(original.dataTable)
      # View(original.dataTable)

      #### TO DO: replace tihs by a parameter sent to iteratetraintest, indicating a subset of features to use ####

      #### Define PC numbers depending on dataset size ####
      stepOver10 <- 10
      pc.numbers <- c(
        2, 3, 4, 5, 6, 7,
        seq(from = 10, to = nrow(original.dataTable) - 1, by = stepOver10),
        nrow(original.dataTable))
      # length(pc.numbers)
#      pc.numbers <- c(2,10,162) ## for quick test

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
                  "\n\tpc.nb:\t\t", pc.nb
          )


          ## Select the first N principal components
          first.pcs <- data.frame(t(dataset$prcomp$x[,1:pc.nb]))
          rownames(first.pcs) <- rownames(t(dataset$prcomp$x[,1:pc.nb]))
          dataset$dataTable <- first.pcs
          # dim(dataset$dataTable)


          ## Define output parameters
          outParam <- outputParameters(
            dataset = dataset,
            classifier = classifier,
            permute = permute,
            createDir = TRUE)

          ## Append suffix with the current number of PCs
          outParam$filePrefix <- paste0(outParam$filePrefix, "_", pc.nb, "_PCs")

          ## Define experiment prefix
          # exp.prefix <-
          #   paste(sep = "_", classifier, dataset$ID, dataset$parameters$feature, dataset$dataType, "nb_of_PCs", pc.nb)
          # if (permute) {
          #   exp.prefix <- paste(sep = "_", exp.prefix,parameters$perm.prefix)
          # }# end if permuted class

          message.with.time("\t", "Experiment prefix: ", outParam$filePrefix)

          #### Run a training/testing experiment ####
          one.result <-  IterateTrainingTesting(
              dataset = dataset,
              classifier = classifier,
              file.prefix = outParam$filePrefix,
              permute = permute
            )
          train.test.results.PCs[[classifier]][[recountID]][[outParam$filePrefix]] <- one.result
          # names(train.test.results.PCs[[classifier]][[recountID]])
        } # end iteration over nb of PCs

      } # end of permutation


      #### Generate boxplot with impact of the number of PCs on the performances of each classifier ####

      ## Reset output parameters
      outParam <- outputParameters(
        dataset = dataset,
        classifier = classifier,
        permute = FALSE,
        createDir = TRUE)
      outParam$filePrefix <- paste0(outParam$filePrefix, "_feature-selection_first-PCs")
      dir.create(path = file.path(outParam$resultDir, "figures"), showWarnings = FALSE, recursive = FALSE)

      ## Generate the error rate boxplot
      ErrorRateBoxPlot(experimentList = train.test.results.PCs[[classifier]][[recountID]],
                       experimentLabels = append(paste(pc.numbers, "PCs"), paste(pc.numbers, "PCs permLabels")),
                       classifier = classifier,
                       horizontal = TRUE,
                       main = paste0(
                         parameters$short_label, " (", parameters$recountID, ") ", parameters$feature, "s",
                         "\n", classifier, "; ", parameters$iterations, " iterations",
                         "\n", "PC selection; ", dataset$dataType),
                       boxplotFile = file.path(
                         outParam$resultDir, "figures",
                         paste(sep = "", outParam$filePrefix, ".pdf")))


      #  train.test.results.all.PCs.per.classifier[[recountID]][[classifier]] <- train.test.results.PCs[[classifier]][[recountID]]
    }  # end of loop over classifiers

    ## Restore original data table
    dataset$dataTable <- original.dataTable
  } # end loop over datasets
} # end loop over recountIDs

## REDO THE ERROR BOXPLOTS
classifier <- classifiers[1]
for (classifier in classifiers) {
  recountID <- selectedRecountIDs[[1]]
  for (recountID in selectedRecountIDs) {
    #  train.test.results.all.PCs.per.classifier[[recountID]] <- list() ## Instantiate an entry per recountID
    parameters <- studyCases[[recountID]]$parameters

    #### Iterate over PCs.to.test ####
    data.type <- "TMM_log2_PC"
    for (data.type in project.parameters$global$PCs.to.test) {
      dataset <- studyCases[[recountID]]$datasetsForTest[[data.type]]

      outParam <- outputParameters(
        dataset = dataset,
        classifier = classifier,
        permute = FALSE,
        createDir = TRUE)
      outParam$filePrefix <- paste0(outParam$filePrefix, "_feature-selection_first-PCs")
      dir.create(path = file.path(outParam$resultDir, "figures"), showWarnings = FALSE, recursive = FALSE)

      ErrorRateBoxPlot(experimentList = train.test.results.PCs[[classifier]][[recountID]],
                       experimentLabels = append(paste(pc.numbers, "PCs"), paste(pc.numbers, "PCs permLabels")),
                       classifier = classifier,
                       horizontal = TRUE,
                       main = paste0(
                         parameters$short_label, " (", parameters$recountID, ") ", parameters$feature, "s",
                         "\n", classifier, "; ", parameters$iterations, " iterations",
                         "\n", "PC selection; ", dataset$dataType),
                       boxplotFile = file.path(
                         outParam$resultDir, "figures",
                         paste(sep = "", outParam$filePrefix, ".pdf")))
    }
  }
}

#### Summarize results for all recountIDs ####
for (classifier in classifiers) {
  allRecountIDsSummary <- list()
  allRecountIDsSummary$perRecountID <- list()
  allRecountIDsSummary$summary <- data.frame()
  # rownames(allRecountIDsSummary$summary) <- names(train.test.results.PCs)
  recountIDs <- names(train.test.results.PCs[[classifier]])
  for (i in 1:length(recountIDs)) {
    recountID <- recountIDs[i]
    oneRecountIDSummary <- SummarizeTrainTestResults(experimentList = train.test.results.PCs[[classifier]][[recountID]])
    allRecountIDsSummary$perRecountID[[recountID]] <- oneRecountIDSummary

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
    paste(collapse = "-", selectedRecountIDs),
    "_", featureType,
    "_feature_selection_first_PCs",
    "_", Sys.Date(), "_results.Rdata"))
message.with.time(
  "Saving results after eval of feature selection with first PCs: ",
  save.result.file)
save(train.test.results.PCs, file = save.result.file)

## Save a memory image that can be re-loaded next time to avoid re-computing all the normalisation and so on.
if (project.parameters$global$save.image) {
  memory.image.file <- file.path(
    project.parameters$global$dir$memoryImages,
    paste0(
      paste(collapse = "-", selectedRecountIDs),
      "_", featureType,
      "_", "feature_selection_first_PCs",
      "_", Sys.Date(), "_memory-image.Rdata"))
  message.with.time("Saving memory image after eval of feature selection with first PCs: ", memory.image.file)
  save.image(file = memory.image.file)
}


message.with.time("Finished script 07_PCA_impact_of_PC_number.R")

