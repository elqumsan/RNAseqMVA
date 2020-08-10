#### Test the impact of k on KNN performances ####

if (length(studyCases) > 1) {
  ## The previous loop over recountIDs is not valid anymore
  stop("Multiple study cases cannot yet be analysed in the same run. \nPlease specify a single value for selected_recount_ids")
}
recountID <- names(studyCases)[1]


message.with.time("Running train/test with all variables to test imapct of knn's parameters for recountID\t", recountID)

## Initialize a list to store all results for the current recountID
train.test.results.all.variables.knn <- list()

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


## Make sure that the classifier variable is KNN
classifier <- "knn"



#### Loop over label permutation ####
permute <- FALSE ## Default for quick test without iterating over all cases
for (permute in project.parameters$global$permute) {

  #### Loop over data types ####
  data.type <- "TMM_log2" ## For test
  for (data.type in project.parameters$global$data.types.to.test) {

    #### Loop over values of k ####
    knn.k <- 3
    for (knn.k in project.parameters$global$knn$k_values) {

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
      outParam <- outputParameters(dataset, classifier, permute, createDir = TRUE, knn.k = knn.k)
      exp.prefix <- outParam$filePrefix
      # exp.prefix <- paste0(recountID,
      #                      "_",
      #                      data.type,
      #                      "_knn_k", knn.k)

      #### Run a training/testing experiment ####
      train.test.results.all.variables.knn[[exp.prefix]] <-
        IterateTrainingTesting(
          dataset,
          classifier = classifier,
          file.prefix = outParam$filePrefix,
          permute = permute)

    } # end iteration over k values

  } # End iterations over datasets

} # End iterations over permutation or not

#### Generate a summary table with the performances ####
summary.table <- data.frame()
permute <- FALSE
for (permute in project.parameters$global$permute) {

  #### Loop over data types ####
  data.type <- "TMM_log2" ## For test
  for (data.type in project.parameters$global$data.types.to.test) {

    #### Loop over values of k ####
    knn.k <- 3
    for (knn.k in project.parameters$global$knn$k_values) {
      dataset <- studyCases[[recountID]]$datasetsForTest[[data.type]]
      dataset$parameters$knn$k <- knn.k

      ## Define output parameters
      outParam <- outputParameters(dataset,
                                   classifier,
                                   permute,
                                   createDir = FALSE,
                                   knn.k = knn.k)
      exp.prefix <- outParam$filePrefix
      current.result <- train.test.results.all.variables.knn[[exp.prefix]]

      current.summary <- data.frame(
        prefix = exp.prefix,
        permute = permute,
        preprocessing = data.type,
        k = knn.k,
        MER = mean(current.result$testing.error.rate)
      )

      if (nrow(summary.table) == 0) {
        summary.table = rbind(summary.table, current.summary)
        #colnames(summary.table) <- c("permute", "normalisation", "k", "MER")
      } else {
        summary.table <- rbind(summary.table, current.summary)
      }
    }
  }
}
# View(summary.table)
# dim(summary.table)

## Compute the summary table
summary.cross.table <- xtabs(MER ~ preprocessing + k + permute, data = summary.table)
# View(summary.cross.table)

cross.table.file <- file.path(
  parameters$dir$tsv,
  paste0(
    recountID,
    "_", featureType,
    "_knn_impact_of_k_result_crosstable.tsv"))
message.with.time("Saving cross-table\t", cross.table.file)
write.table(x = as.matrix(summary.cross.table[,,1]),
            file = cross.table.file,
            quote = FALSE, sep = "\t",
            row.names = TRUE, col.names = NA)

cross.table.file.perm <- file.path(
  parameters$dir$tsv,
  paste0(
    recountID,
    "_", featureType,
    "_knn_impact_of_k_result_crosstable_permLabels.tsv"))
message.with.time("Saving cross-table with permutation test result\t", cross.table.file)
write.table(x = as.matrix(summary.cross.table[,,2]),
            file = cross.table.file.perm,
            quote = FALSE, sep = "\t",
            row.names = TRUE, col.names = NA)

# library(knitr)
# kable(summary.cross.table[,,1])

#### Plotting the Misclassification Error rate using all diverse data type all variables with KNN classifier? ####
outParam <- outputParameters(
  dataset,
  classifier = "knn impact of k",
  permute = FALSE, createDir = TRUE)
dir.create(path = file.path(outParam$resultDir, "figures"), showWarnings = FALSE, recursive = FALSE)


## Suppress recount ID and feature type from the labels, since they are always the same for one experiment
experimentLabels <- names(train.test.results.all.variables.knn)
experimentLabels <- sub(pattern = paste0(parameters$recountID, "_",
                                         parameters$feature, "_knn_"),
                        replacement = "",
                        x = experimentLabels)
## Replace underscore by space in the labels
experimentLabels <- sub(pattern = "_", replacement = " ", x = experimentLabels)

## TEMPORARY FIX FOR THE ORDER OF THE EXPERIMENTS
k <- sub(pattern = " .*", replacement = "", x = experimentLabels, perl = TRUE)
dataType <- sub(pattern = "[^ ]+ ", replacement = "", x = experimentLabels, perl = TRUE)
#  nonPermuted <- grep(pattern = "permLabels", x = experimentLabels, invert = TRUE)
#  experimentOrder <- order(dataType[nonPermuted])
experimentOrder <- order(dataType)


ErrorRateBoxPlot(experimentList = train.test.results.all.variables.knn[experimentOrder],
                 classifier = classifier,
                 horizontal = TRUE,
                 experimentLabels = experimentLabels[experimentOrder],
                 # boxplotFile = file.path(
                 #   outParam$resultDir, "figures",
                 #   paste(sep = "", outParam$filePrefix, ".pdf")),
                 boxplotFile = file.path(
                   outParam$resultDir, "figures",
                   paste0(parameters$recountID,
                          "_", parameters$feature,
                          "_knn_impact_of_k",
                          ".pdf")),
                 main = paste(sep = "",
                              parameters$short_label, " (", parameters$recountID, ") ",
                              parameters$feature,
                              "\nImpact of k on KNN; ",
                              project.parameters$global$iterations, " iterations")
)


## Save the results in a separate object, that can be reloaded later
## Define the path to the memory image for this test (compare classifier whenn they use all variables as features)
save.result.file <- file.path(
  project.parameters$global$dir$memoryImages,
  paste0(
    paste(collapse = "-", selectedRecountIDs),
    "_", featureType,
    "_knn_impact_of_k_result.Rdata"))
dir.create(project.parameters$global$dir$memoryImages, showWarnings = FALSE, recursive = TRUE)
message.with.time("Saving results  after eval of impact of k on knn: ", save.result.file)
save(train.test.results.all.variables.knn, file = save.result.file)

## Save a memory image that can be re-loaded next time to avoid re-computing all the normalisation and so on.
if (project.parameters$global$save.image) {
  memory.image.file <- file.path(
    project.parameters$global$dir$memoryImages,
    paste0(
      paste(collapse = "-", selectedRecountIDs),
      "_", featureType,
      "_", "knn_impact_of_k_memory-image.Rdata"))
  message.with.time("Saving memory image after eval of impact of k on knn: ", memory.image.file)
  save.image(file = memory.image.file)
}
memory.image.file

message.with.time("Finished script 13_knn_impact_of_k.R")
