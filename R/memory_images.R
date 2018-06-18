
#' @title Save a memory image in order to be able to reload all results without having to recompute everything
#' @author Jacques van Helden anf Mustafa AbuElQumsan
#' @description if you want to saving time without of made recomputation thoughout loading last computation process.
#' @param memory.dir = parameters$dir$memoryImages it is the path the saving the image file.
#' @param recountID = parameters$recountID Recount ID number.
#' @param prefix = "RNAseqMVA_image_" such prefix in order to saving the file in this prefix.
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
#' @description if you want to loading abeforehand image file which is prepared from the saveMemoryImage ()
#' @param memory.dir = parameters$dir$memoryImages it is the path the saving the image file.
#' @param recountID = parameters$recountID Recount ID number.
#' @param prefix = "RNAseqMVA_image_" such prefix in order to saving the file in this prefix.
#' @export
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

