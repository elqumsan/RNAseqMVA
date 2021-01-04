#### Reload SVM results from the different study cases in order to compare performances ####

#### Initialise variable to store results ####
kernel.impact.results <- list()
# MERmean <- list()
# MERsd <- list()
# MERiter <- list()
MERstat <- list()
outfiles <- list()
iterations <- NULL

#### Parameter values ####

## Output directories
dir <- list()
dir[["figures"]] <- project.parameters$global$dir$figures
dir[["tables"]] <- project.parameters$global$dir$TSV

## Graphical parameters
par.ori <- par(no.readonly = TRUE)

## Feature types
featureTypes <- c("transcript", "gene")

## Study case IDs
# recountIDs <- c(
#   "SRP042620",  # Breast cancer
#   "SRP062966",  # Lupus (sc)
#   "SRP035988",  # Psoriasis
#   "SRP056295",  # Acute myeloid leukemia
#   "SRP057196",  # Brain cells (sc)
#   "SRP066834"   # Cerebral organoids and fetal neocortex (sc)
# )

studyCaseLabels <- c(
  "SRP042620" = "Breast cancer",
  "SRP056295" = "Acute myeloid leukemia",
  "SRP035988" = "Psoriasis",
  "SRP057196" = "Brain cells (sc)",
  "SRP066834" = "Cerebral organoids and fetal neocortex (sc)",
  "SRP062966" = "Lupus (sc)"
)

recountIDs <- names(studyCaseLabels)

## Kernels
kernels <- c(
  "linear",
  "polynomial",
  "sigmoid",
  "radial"
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

## Make a model for the stat matrices (will be used for different stats: mean, median, sd, ...)
statDataFrame <- data.frame(
  matrix(nrow = length(normMethods),
         ncol = length(kernels) * length(permSuffixes)))
rownames(statDataFrame) <- normMethods
colnames(statDataFrame) <-
      apply(expand.grid(kernels, permSuffixes), 1, paste, collapse="")

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
    for (stat in c("mean", "median", "sd", "min", "Q1", "median", "Q3", "max")) {
      MERstat[[featureType]][[recountID]][[stat]] <- statDataFrame
    }
    # MERstat[[featureType]][[recountID]] <-
    #   data.frame(
    #     matrix(nrow = length(normMethods),
    #            ncol = length(kernels) * length(permSuffixes)))
    #
    # ## Misclassification Error Rate (MER)
    # MERmean[[featureType]][[recountID]] <-
    #   data.frame(matrix(nrow = length(normMethods),
    #                     ncol = length(kernels) * length(permSuffixes)))
    # rownames(MERmean[[featureType]][[recountID]]) <- normMethods
    # colnames(MERmean[[featureType]][[recountID]]) <-
    #   apply(expand.grid(kernels, permSuffixes), 1, paste, collapse="")
    #
    # ## Standard deviation of the MER
    # MERsd[[featureType]][[recountID]] <-
    #   data.frame(matrix(nrow = length(normMethods),
    #                     ncol = length(kernels) * length(permSuffixes)))
    # rownames(MERsd[[featureType]][[recountID]]) <- normMethods
    # colnames(MERsd[[featureType]][[recountID]]) <-
    #   apply(expand.grid(kernels, permSuffixes), 1, paste, collapse="")
    #
    # ## Number of iterations for each test
    # MERiter[[featureType]][[recountID]] <-
    #   data.frame(matrix(nrow = length(normMethods),
    #                     ncol = length(kernels) * length(permSuffixes)))
    # rownames(MERiter[[featureType]][[recountID]]) <- normMethods
    # colnames(MERiter[[featureType]][[recountID]]) <-
    #   apply(expand.grid(kernels, permSuffixes), 1, paste, collapse="")
    #
    # ## Median
    # MERmedian[[featureType]][[recountID]] <-
    #   data.frame(matrix(nrow = length(normMethods),
    #                     ncol = length(kernels) * length(permSuffixes)))
    # rownames(MERmedian[[featureType]][[recountID]]) <- normMethods
    # colnames(MERmedian[[featureType]][[recountID]]) <-
    #   apply(expand.grid(kernels, permSuffixes), 1, paste, collapse="")

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
            MERstat[[featureType]][[recountID]][["mean"]][normMethod, paste0(kernel, permSuffix)] <-
              mean(result$testing.error.rate)
            MERstat[[featureType]][[recountID]][["sd"]][normMethod, paste0(kernel, permSuffix)] <-
              sd(result$testing.error.rate)
            MERstat[[featureType]][[recountID]][["min"]][normMethod, paste0(kernel, permSuffix)] <-
              min(result$testing.error.rate)
            MERstat[[featureType]][[recountID]][["Q1"]][normMethod, paste0(kernel, permSuffix)] <-
              quantile(result$testing.error.rate, prob = 0.25)
            MERstat[[featureType]][[recountID]][["median"]][normMethod, paste0(kernel, permSuffix)] <-
              median(result$testing.error.rate)
            MERstat[[featureType]][[recountID]][["Q3"]][normMethod, paste0(kernel, permSuffix)] <-
              quantile(result$testing.error.rate, prob = 0.75)
            MERstat[[featureType]][[recountID]][["max"]][normMethod, paste0(kernel, permSuffix)] <-
              max(result$testing.error.rate)
          } else {
            warning("Missing result\t", resultName)
          }
        }
      }
    }
  }
}



### Barplots ####
featureType <- "gene"
for (featureType in featureTypes) {
  outfiles[[paste0(featureType, "_barplots")]] <- file.path(
    dir[["figures"]],
    paste0("svm_kernel-impact_barplots_", featureType, ".pdf"))

  message("Exporting barplots in file ",  outfiles[[paste0(featureType, "_barplots")]])
  pdf(file = outfiles[[paste0(featureType, "_barplots")]],
      width = 8,
      height = 8)
  par(mfrow  = c(3, 2))
  par(mar = c(6, 5, 4, 1))

  for (recountID in recountIDs) {
    ## Collect stat matrices
    MERmean <- t(as.matrix(MERstat[[featureType]][[recountID]][["mean"]][, 1:length(kernels)]))
    MERmin <- t(as.matrix(MERstat[[featureType]][[recountID]][["min"]][, 1:length(kernels)]))
    MERmax <- t(as.matrix(MERstat[[featureType]][[recountID]][["max"]][, 1:length(kernels)]))
    MERq1 <- t(as.matrix(MERstat[[featureType]][[recountID]][["Q1"]][, 1:length(kernels)]))
    MERq3 <- t(as.matrix(MERstat[[featureType]][[recountID]][["Q3"]][, 1:length(kernels)]))

    z <- barplot(MERmean,
            beside = TRUE, horiz =FALSE,
            main = paste0(studyCaseLabels[recountID], "\n(", recountID, "); ", featureType),
            las = 2,
            xlim = c(0, length(normMethods) * length(kernels) * 1.65),
            ylim = c(0,0.8),
            cex.names = 0.8,
            border = "#FFFFFF",
            ylab = "MER",
            legend.text = kernels,
            args.legend = list(cex = 0.7))
    arrows(x0 = z,
           y0 = MERmin,
           x1 = z,
           y1 = MERmax,
           angle=90, code=3, length=0.0)
    rect(xleft = as.vector(z - 0.2),
         ybottom = as.vector(MERq1),
         xright = as.vector(z + 0.2),
         ytop = as.vector(MERq3),
         col = "white",
         border = "black")
  }
  par(mfrow = c(1,1))
  par(par.ori)
  dev.off()
}


library(knitr)
kable(t(as.data.frame(outfiles)), col.names = "File")

