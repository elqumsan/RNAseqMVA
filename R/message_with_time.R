#' @title print a message preceded by currrent tume
#' @author Jacques van Helden
#' @param ... all parameters are passed to the message() function
#' @export
message.with.time <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d_%H%M%S"), "\t", ...)
}


