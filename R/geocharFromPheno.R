#' @title Extract the conditions from the "characteristics" column of the coldata.
#' @author Mustafa ABUELQUMSAN and Jacques van Helden
#' @description This is a bit tricky: we have to parse a string describing several attributes.
#' @param  runPheno it is the Run pheno table which may be containt duplicate in the rus for the same sample, and it give us the actual length for a pheno table
#' @import recount
#' @import SummarizedExperiment
#'
#' @export
geocharFromPheno <- function(runPheno) {
  message("Extracting geo characteristics")
  geochar.list <- lapply(split(
    runPheno,
    seq(from=1, to=nrow(runPheno))),
    geo_characteristics)
  # class(geochar.list)
  # View(geochar.list)
  # names(geochar.list)
  # head(geochar.list)

  ## collect column names
  geo_colnames <- vector()
  for (i in 1:length(geochar.list)) {
    geo_colnames <- union(geo_colnames, names((geochar.list[i])[[1]]))
  }

  ## Convert the list to a data.frame
  geochar.frame <- data.frame(matrix(nrow = length(geochar.list), ncol=length( geo_colnames))) ## Instantiate the data frame with the first sample
  names(geochar.frame) <- geo_colnames
  rownames(geochar.frame) <- names(geochar.list)
  i <- 27
  for (i in 1:length(geochar.list)) {
    sample.info <- (geochar.list[i])[[1]]
    geochar.frame[i, colnames(sample.info)] <- sample.info
  }

  ## Replace column header "cells" by "cell.line"
  ## This code comes from Leonardo, but we don't think we really need them for our purposes
  ## The rbind does not work for some samples (e.g. SRP006574) because the number of columns of geochar varies from sample to sample.
  # geochar <- do.call(rbind,
  #                    lapply(geochar.list,
  #                                  function(x) {
  #                                    if ('cells' %in% colnames(x)) {
  #                                      colnames(x)[colnames(x) == 'cells'] <- 'cell.line' ## Replace column header "cells" by "cell.line"
  #                                      return(x)
  #                                    } else {
  #                                      return(x)
  #                                    }
  #                                  }))

  return(geochar.frame)
}
