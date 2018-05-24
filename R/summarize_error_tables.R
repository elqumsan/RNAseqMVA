#' @title generate a summary table from a list of error table files.
#' @description read a list of error tables and summarize their results in
#' a single summary table with 1 row per input file and one column per statistics
#' of the original error tables. The summary statistics is the mean, except for
#' the iteration, where the max value of the error table is taken.
#' @author Jacques van Helden and Mustafa AbuElQumsan
#' @param files a list of files with the error rates for one classification analysis.
#' @param stopIfMissing=TRUE if TRUE, issue an error if any of the input files does not exist
#' @param labels=files labels associated to each file (used for displays)
#' The files should be formatted as tab-separated values.
#' Each row corresponds to one iteration.
#'
#' @import recount
#'
#' @export
SummarizeErrorTable <- function (files,
                                 stopIfMissing = TRUE,
                                 labels = files) {

  ## Check if there are missing files
  missing.files <- files[!file.exists(files)]
  if (length(missing.files) > 0) {
    message("Missing files\n\t", paste(collapse="\n\t", missing.files))
    if (stopIfMissing) {
      stop("ErrorRateComparativePlot() stopped because ", length(missing.files), " files are missing")
    }
  }


  ## Check labels
  if (length(labels) != length(files)) {
    stop("ErrorRateComparativePlot() arguments 'files' and 'labels' must have the same length")
  }


  #### Read the error tables ####
  checked.files <- files[file.exists(files)]
  message("Reading ", length(checked.files), " error tables. ")
#  errors <- list()
  summary.table <- data.frame()
  for (file in files) {
    if (file.exists(file)) {

      current.errors <- read.delim(file = file, header = TRUE, sep = "\t")
      names(current.errors)[1] <- "iterations"

      summary.row <- apply (current.errors, 2, mean)
      summary.row["iterations"] <- max(current.errors$iteration)

#       if (nrow(summary.table) == 0) {
#         summary.table <- t(data.frame(summary.row))
# #        row.names(summary.table) <- file
#       } else {
#         summary.table[file,] <- summary.row
#       }
      summary.table <- rbind(summary.table, summary.row)
      names(summary.table) <- names(current.errors)

#      errors[[file]] <- current.errors
      rm(current.errors)

    } else {
      summary.table <- rbind(summary.table, rep(NA, length.out = ncol(summary.table)))

    }
  }

  summary.table$label <- files
  summary.table$file <- files
  # View(summary.table)
  return(summary.table)
}
