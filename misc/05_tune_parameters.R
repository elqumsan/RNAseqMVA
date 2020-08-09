#### Use e1071::tune() wrapper to find optimal parameters for the different classifiers ####

## Note: e1071::tune() functions are versy convenient to analyse user-specified combinations of
## parameters but they run all the tests sequentially, which may take a huge time for big datasets
## (as the ones we analyse) and if the number of parameter values is high.
##
## Moreover, the tuning tests each parameter value only once, but we noticed that for some parameters
## there are strong fluctuations between successive runs (e.g. the ntrees and mtry parameters of randomForest).
## For this reason, we developed a specific procedure to run evaluate the impact of parameters with a
## user-specified number of iterations (which enables us to evaluate the fluctuatins of performances
## between independent runs).
##
## We however keep this code for the sake of comparison.

## Select a data type for the tuning
dataType <- "TMM_log2"


## Define the path to the memory image for this test (compare classifier whenn they use all variables as features)
dir.create(project.parameters$global$dir$memoryImages, showWarnings = FALSE, recursive = TRUE)
featureType <- project.parameters$global$feature
save.result.file <- file.path(
  project.parameters$global$dir$memoryImages,
  paste0(
    studyCase$ID,
    "_", featureType,
    "_", "tune-classifiers",
    "_", dataType,
    "_results.Rdata"))
# load(save.result.file)

#### Run the analysis ####
tuningResult <- TuneClassifiers(
  studyCase = studyCase,
  dataType = dataType,
  tuneSVM = TRUE,
  tuneRandomForest = TRUE,
  tuneKNN = TRUE,
  plotResults = TRUE)


#### Export plots ####

## Define directory and file prefix for the plots
resultDir <- file.path(studyCase$parameters$dir$results, "tuning")
dir.create(resultDir, showWarnings = FALSE, recursive = TRUE)
filePrefix <- file.path(resultDir,
                        paste0(
                          studyCase$ID,
                          "_", studyCase$parameters$feature,
                          "_", dataType,
                          "_tuning"))

plotRFtuning(tuningResult, ylim = c(0,1), pdfFile = paste0(filePrefix, "_RF.pdf"))
plotKNNtuning(tuningResult, ylim = c(0,1), pdfFile = paste0(filePrefix, "_KNN.pdf"))


#### Save the results in a separate object, that can be reloaded later ####
message.with.time(
  "Saving results after parameter tuning: ",
  save.result.file)
save(tuningResult, file = save.result.file)

## Ending message
message.with.time("Finished script misc/05_tune_parameters.R")

