################################################################
## Run training-testing iterations with increasing number of principal components.

##### define the file to store memory Image for " the Number of PCs" test #####
image.dir <- file.path( parameters$dir$memoryImages, parameters$recountID)
dir.create(image.dir, showWarnings = FALSE, recursive = TRUE)

image.file <- file.path(
  project.parameters$global$dir$memoryImages,
  paste(sep = "", "classif_eval_nb-of-PCs_",
        paste(collapse = "-", selectedRecountIDs),
        "_", Sys.Date(), ".Rdata"))


#image.file <- file.path(image.dir, paste(sep = "", "train_test_no._of_PCs_", parameters$recountID , ".Rdata"))


## For debug: reset the parameteres for all the study cases
## This is used to re-run the analyses on each study case after having changed some parameters in the yaml-specific configuration file
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
recountID <- selectedRecountIDs[[1]]
for (recountID in selectedRecountIDs) {
  message("Impact of PCs\t", recountID)


  #  train.test.results.all.PCs.per.classifier[[recountID]] <- list() ## Instantiate an entry per recountID
  parameters <- studyCases[[recountID]]$parameters
  #studyCases$PCsVar <- PCsWithTrainTestSets(studyCases$filtered)

  classifier <- "svm" ## For testing
  for (classifier in project.parameters$global$classifiers) {
    ## List to store all results
    train.test.results.PCs[[recountID]] <- list()


    #### Iterate over PCs.to.test ####
    data.type <- "TMM_log2_PC"
    for (data.type in project.parameters$global$PCs.to.test) {

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

          ## define experiment prefix
          currentOutputParameters <- outputParameters(dataset = dataset, classifier = classifier, permute = permute, createDir = TRUE)
          exp.prefix <-
            paste(sep = "_", classifier, dataset$ID  , dataset$dataType, "nb_of_PCs", pc.nb)
          if (permute) {
            exp.prefix <- paste(sep = "_", exp.prefix,parameters$perm.prefix)
          }# end if permuted class

          message.with.time("\t", "Experiment prefix: ", exp.prefix)

          train.test.results.PCs[[recountID]][[exp.prefix]] <-
            IterateTrainingTesting(
              dataset,
              classifier = classifier,
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
      ErrorRateBoxPlot(experimentList = train.test.results.PCs[[recountID]],
                       classifier = classifier,
                       horizontal = TRUE,
                       main = paste0(
                         "PC selection; ", parameters$recountID, "; ", dataset$dataType, "\n",
                         classifier, "; ", parameters$iterations , "iterations,"),
                       boxplotFile = file.path(
                         outParam$resultDir, "figures",
                         paste(sep = "", outParam$filePrefix, ".pdf"))
      )


      #  train.test.results.all.PCs.per.classifier[[recountID]][[classifier]] <- train.test.results.PCs[[recountID]]
    }  # end of loop over classifiers

    ## Restore original data table
    dataset$dataTable <- original.dataTable
  } # end loop over datasets
} # end loop over recountIDs

