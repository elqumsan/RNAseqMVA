#### Reload SVM results from the different study cases in order to compare performances ####

## Initialise variable to store results
kernel.impact.results <- list()

## Parameter values
featureTypes <- c("gene", "transcript")
recountIDs <- c(
  "SRP042620",
  "SRP062966",
  "SRP042620",
  "SRP035988",
  "SRP056295",
  "SRP057196",
  "SRP066834"
)

## Default values for testing
featureType <- "transcript"
recountID <- "SRP042620"

for (featureType in featureTypes) {
  for (recountID in recountIDs) {

    ## Define the path to the memory image for this test (compare classifier whenn they use all variables as features)
    save.result.file <- file.path(
      project.parameters$global$dir$memoryImages,
      paste0(
        recountID,
        "_", featureType,
        "_svm_impact_of_kernel",
        "_results.Rdata"))

    if (file.exists(save.result.file)) {
      message("Reloading results for study case\t", recountID, "\tfeature type\t", featureType)
      load(save.result.file)
      kernel.impact.results[[featureType]][[recountID]] <- train.test.results.all.variables.per.svm[[recountID]]
    } else {
      message("Missing result file\t", save.result.file)
    }
  }
}

## Check the names of the loaded study cases
for (featureType in featureTypes) {
  message("Loaded study cases for feature type ", featureType)
  message("\t", paste(collapse = ", ", names(kernel.impact.results[[featureType]])))
}



