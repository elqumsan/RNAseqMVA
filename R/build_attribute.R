#' @title build attributes for an object depending on its class
#' @export
buildAttributes <- function(self) {
  message("\tBuilding attributes for object of class ", paste(collapse=", ", class(self)))
  UseMethod("buildAttributes", self)
#  message("trainIndices length : ", length(self$trainTestProperties$trainIndices))
#  message("\t\treturning from buildAttributes()")
  return(self) ## NOT SURE THIS LINE IS REQUIRED. TO BE CHECKED.
}

#' @title default method to build attributes for an object.
#' @description Just send message with object classes. The class-specific builders should have been be called before.
#' @export
buildAttributes.default <- function(self) {
  message("\tFinished building attributes for object of class ", paste(collapse=",", class(self)))
#  print (self$trainTestProperties)
#  message("trainIndices length : ", length(self$trainTestProperties$trainIndices))
#  message("\t\treturning from buildAttributes.default()")
  return(self)
}


