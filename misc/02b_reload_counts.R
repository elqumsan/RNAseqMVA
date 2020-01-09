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

