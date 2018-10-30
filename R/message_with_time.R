#' @title print a message preceded by currrent time
#' @author Jacques van Helden and Mustafa AbuElQumsan
#' @description Print a message prepended with the current time in reverse order (YYYY-MM-DD) in order to match alphabetical order with chronological order.
#' @param ... all parameters are passed to the message() function
#'
#' @return report message with time for each proccess.
#'
#' @export
message.with.time <- function (...) {
  message(format(Sys.time(), "%Y-%m-%d_%H%M%S"), "\t", ...)
}


