#### Reload SVM results from the different study cases in order to compare performances ####

#### Initialise variable to store results ####
kernel.impact.results <- list()
testingMERmean <- list()
testingMERsd <- list()

#### Parameter values ####
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


## Kernels
kernels <- c(
  "radial",
  "sigmoid",
  "polynomial",
  "linear"
)

## Normalisation methods
normMethods <- c(
  "filtered_counts",
  "q0.75",
  "TMM",
  "RLE",
  "q0.75_log2",
  "TMM_log2",
  "RLE_log2",
  "q0.75_log2_PC",
  "TMM_log2_PC",
  "RLE_log2_PC"
)

## Permutation test suffixes
permSuffixes <- c("", "_permLabels")

#### Default values for debugging and testing ####
featureType <- "gene"
recountID <- "SRP042620"
normMethod <- "TMM_log2"
kernel <- "linear"
permSuffix <- ""



#### Reload the results ####
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

#### Compare performances between kernels for each study case ####
for (featureType in featureTypes) {
  for (recountID in recountIDs) {

    ## Initialise the comparison tables
    testingMERmean[[featureType]][[recountID]] <-
      data.frame(matrix(nrow = length(normMethods),
                        ncol = length(kernels) * length(permSuffixes)))
    rownames(testingMERmean[[featureType]][[recountID]]) <- normMethods
    colnames(testingMERmean[[featureType]][[recountID]]) <-
      apply(expand.grid(kernels, permSuffixes), 1, paste, collapse="")

    testingMERsd[[featureType]][[recountID]] <-
      data.frame(matrix(nrow = length(normMethods),
                        ncol = length(kernels) * length(permSuffixes)))
    rownames(testingMERsd[[featureType]][[recountID]]) <- normMethods
    colnames(testingMERsd[[featureType]][[recountID]]) <-
      apply(expand.grid(kernels, permSuffixes), 1, paste, collapse="")

    for (permSuffix in permSuffixes) {
      for (normMethod in normMethods) {
        for (kernel in kernels) {
          resultName <- paste(sep="_",
                              recountID,
                              featureType,
                              "svm",
                              kernel,
                              normMethod)
          resultName <- paste0(resultName, permSuffix)
          if (resultName %in% names(kernel.impact.results[[featureType]][[recountID]])) {
            result <- kernel.impact.results[[featureType]][[recountID]][[resultName]]
            # names(result)
            # View(result)
            ## Compute mean misclassification error rate
            testingMERmean[[featureType]][[recountID]][normMethod, paste0(kernel, permSuffix)] <-
              mean(result$testing.error.rate)
            testingMERsd[[featureType]][[recountID]][normMethod, paste0(kernel, permSuffix)] <-
              sd(result$testing.error.rate)
          } else {
            warning("Missing result\t", resultName)
          }
        }
      }
    }
  }
}

testingMERmean[[featureType]][[recountID]]


testingMERsd[[featureType]][[recountID]]


