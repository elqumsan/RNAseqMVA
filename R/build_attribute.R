#' @title build attributes for an object depending on its class
#' @export
buildAttributes <- function(self) {
  message("\tBuilding attributes for object of class ", paste(collapse=", ", class(self)))
  UseMethod("buildAttributes", self)
  message("returning from buildAttributes")
  return(self)
}

#' @title default method to build attributes for an object.
#' @description Just send message with object classes. The class-specific builders should have been be called before.
#' @export
buildAttributes.default <- function(self) {
  message("\tFinished building attributes for object of class ", paste(collapse=",", class(self)))
#  print (self$trainTestProperties)
  message("returning from buildAttributes.default")
  return(self)
}

