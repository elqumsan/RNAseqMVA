#' @title print a message preceded by currrent tume
#' @author Jacques van Helden and Mustafa AbuElQumsan
#' @description such function in order to be able to compute the time has been exploit through the computation process.
#' @param ... all parameters are passed to the message() function
#' @return report message with time for each proccess.
#' @export
message.with.time <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d_%H%M%S"), "\t", ...)
}


