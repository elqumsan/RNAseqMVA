# ==============================================================================
# Parameter Tuning for Machine Learning Classifiers
# ==============================================================================
# This script uses e1071::tune() wrapper to find optimal parameters for
# different classifiers (SVM, Random Forest, KNN).
#
# Note: e1071::tune() functions are very convenient for analyzing user-specified
# combinations of parameters, but they run all tests sequentially, which may
# take significant time for large datasets and high numbers of parameter values.
#
# Additionally, the tuning tests each parameter value only once, but we have
# observed strong fluctuations between successive runs for some parameters
# (e.g., ntrees and mtry parameters of randomForest). For this reason, we
# developed a specific procedure to evaluate parameter impact with a
# user-specified number of iterations to assess performance fluctuations
# between independent runs.
#
# This code is kept for comparison purposes.
# ==============================================================================

# Configuration ----------------------------------------------------------------

# Select data type for tuning
dataType <- "TMM_log2"

# Define feature type
featureType <- project.parameters$global$feature

# Setup directories and file paths --------------------------------------------

# Create memory images directory if it doesn't exist
dir.create(project.parameters$global$dir$memoryImages,
           showWarnings = FALSE,
           recursive = TRUE)

# Define path for saving results
save.result.file <- file.path(
  project.parameters$global$dir$memoryImages,
  paste0(
    studyCase$ID,
    "_", featureType,
    "_", "tune-classifiers",
    "_", dataType,
    "_results.Rdata"
  )
)

# Uncomment to load existing results
# load(save.result.file)

# Run parameter tuning analysis -----------------------------------------------

message.with.time("Starting parameter tuning analysis...")

tuningResult <- TuneClassifiers(
  studyCase = studyCase,
  dataType = dataType,
  tuneSVM = TRUE,
  tuneRandomForest = TRUE,
  tuneKNN = TRUE,
  plotResults = TRUE
)

# Export visualization plots --------------------------------------------------

# Define output directory and file prefix for plots
resultDir <- file.path(studyCase$parameters$dir$results, "tuning")
dir.create(resultDir, showWarnings = FALSE, recursive = TRUE)

filePrefix <- file.path(
  resultDir,
  paste0(
    studyCase$ID,
    "_", studyCase$parameters$feature,
    "_", dataType,
    "_tuning"
  )
)

# Generate tuning plots with consistent formatting
plotRFtuning(tuningResult,
             ylim = c(0, 1),
             pdfFile = paste0(filePrefix, "_RF.pdf"))

plotKNNtuning(tuningResult,
              ylim = c(0, 1),
              pdfFile = paste0(filePrefix, "_KNN.pdf"))

# Save results -----------------------------------------------------------------

message.with.time(
  "Saving results after parameter tuning: ",
  save.result.file
)

save(tuningResult, file = save.result.file)

# Completion message
message.with.time("Finished script misc/05_tune_parameters.R")
