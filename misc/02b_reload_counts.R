#################### Reload previously saved and normalized counts ####################

## Define the path of the memory image file
featureType <- project.parameters$global$feature
studyCases.mem.image <- file.path(
  project.parameters$global$dir$memoryImages,
  paste0(
    paste(collapse = "-", selectedRecountIDs),
    "_", featureType,
    "_loaded_studyCases.Rdata"))

## Reload previously stored memory image
message.with.time("Reloading study cases from previously stored memory image",
        "\n\t", studyCases.mem.image)
system.time(load(studyCases.mem.image))

## Indicate that this script has finished running
message.with.time("Finished running 02b_reload_counts.R")

## Print out some statistics about the data sizes
message("Data set characteristics")

## Extract parameters for the current study case
studyCase <- studyCases[[recountID]]
# View(studyCase)
message("Computing filtering impact for study case ", recountID)
message("Recount ID: ", studyCase$rawData$countsPerSample$parameters$recountID)
message("Name: ", studyCase$rawData$countsPerSample$parameters$short_label)

datasetParameters <- data.frame(
  beforeFiltering = c(
    nbClasses = studyCase$rawData$countsPerSample$nbClasses,
    nbSamples = studyCase$rawData$countsPerSample$nbSamples,
    nbGenes = studyCase$rawData$countsPerSample$nbGenes
  ),
  afterFiltering = c(
    nbClasses = studyCase$datasetsForTest$filtered$nbClasses,
    nbSamples = studyCase$datasetsForTest$filtered$nbSamples,
    nbGenes = studyCase$datasetsForTest$filtered$nbGenes
  )
)

discarded_feature = c(
  nbGenes = studyCase$rawData$countsPerSample$nbGenes,
  kept_after_NZVF = length(studyCase$datasetsForTest$filtered$keptGenes),
  discarded_by_NZVF = studyCase$rawData$countsPerSample$nbGenes - length(studyCase$datasetsForTest$filtered$keptGenes)
  )

datasetParameters$discarded <- datasetParameters$beforeFiltering - datasetParameters$afterFiltering
print(datasetParameters)
print(discarded_feature )
