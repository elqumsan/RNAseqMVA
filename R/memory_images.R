
#' @title Save a memory image in order to be able to reload all results without having to recompute everything
#' @author Jacques van Helden anf Mustafa AbuElQumsan
#' @export
saveMemoryImage <- function(memory.dir = parameters$dir$memoryImages,
                           recountID = parameters$recountID,
                           prefix = "RNAseqMVA_image_") {
  message.with.time("Saving memory image")
  image.dir <- file.path (memory.dir, recountID)
  dir.create(image.dir, showWarnings = FALSE, recursive = TRUE)
  image.file <- file.path(image.dir,  paste(sep="", prefix, recountID,".Rdata"))
  message("\tMemory image file: ", image.file)
  save.image(file = image.file)
}

#' @title Load a previous analysis previously stored as memoy image
#' @author Jacques van Helden anf Mustafa AbuElQumsan
#' @export
loadMemoryImage <- function(memory.dir = parameters$dir$memoryImages,
                            recountID = parameters$recountID,
                            prefix = "RNAseqMVA_image_") {
  message.with.time("Loading memory image")
  image.dir <- file.path (memory.dir, recountID)
  image.file <- file.path(image.dir,  paste(sep="", prefix, recountID,".Rdata"))
  if (file.exists(image.file)) {
    message("\tMemory image file: ", image.file)
    load(file = image.file)
  } else {
    stop("Cannot reload memory image file ", image.file)
  }
}

