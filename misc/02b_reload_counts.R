#################### Reload previously loaded and normalized counts ####################

## Reload previously stored memory image
if (!is.null(project.parameters$global$reload.date)) {
  image.date <- project.parameters$global$reload.date
} else {
  image.date <- Sys.Date()
}
featureType <- project.parameters$global$feature
studyCases.mem.image <- file.path(
  project.parameters$global$dir$memoryImages,
  paste(sep = "", "loaded_studyCases_",
        paste(collapse = "-", selectedRecountIDs),
        "_", featureType,
        "_", image.date, ".Rdata"))

message("Reloading study cases from previously stored memory image",
        "\n\t", studyCases.mem.image)
load(studyCases.mem.image)



## Indicate that this script has finished running
message.with.time("Finished running 02b_reload_counts.R")

