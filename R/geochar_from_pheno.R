#' @title Extract the conditions from the "characteristics" column of Recount pheno Table, and store them in separate columns.
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @description This is a bit tricky: we have to parse a string describing several attributes.
#' @param  runPheno run pheno table produced by colData(), which contains the pheno information for each run (there might be several runs per )
#' @return phenoTable with data.frame properities in order to be much earier in handle it smoothly.
#'
#' @export
geocharFromPheno <- function(runPheno) {

  LoadRequiredBioconductorPackages(c("recount", "SummarizedExperiment"))

  message("\tExtracting geo characteristicsfrom Recount pheno table")

  ### TEMPORARY PATCH FOR BUG IN THE RECOUNT
  ###  PHENO TABLE OF SOME TRANSCRIPT DATASET
  if (class(runPheno$characteristics) == "character") {
    geo.characteristics <-
      as.data.frame(t(as.data.frame(
      strsplit(runPheno$characteristics,
               split = "; *",
               perl = TRUE))))
    rownames(geo.characteristics) <- NULL

    ## Parse column names
    c <- 1
    for (c in 1:ncol(geo.characteristics)) {
      first.value <- as.vector(unlist((geo.characteristics[1, 1])))
      first.value.split <- unlist(strsplit(first.value, split = ":. *", perl = TRUE))
      if (length(first.value.split) == 2) {
        field.name <- first.value.split[1]
        names(geo.characteristics)[c] <-  field.name
        geo.characteristics[c] <- sub(
          x = geo.characteristics[,c],
          pattern = paste0(field.name, ": *"),
          replacement = "", perl = TRUE)
      }
    }
  } else {
    geo.characteristics <- recount::geo_characteristics(runPheno)
  }
  # View(geo.characteristics)
  # class(geo.characteristics)
  # names(geo.characteristics)
  # dim(geo.characteristics)
  # table(geo.characteristics)
  message("\t\t", nrow(runPheno), " rows (runs). ")
  message("\t\tgeo_characteristics fields found\t",
          paste(collapse = ", ", colnames(geo.characteristics)))

  return(geo.characteristics)

  # View(geo.characteristics)
  # geochar.list <- lapply(split(
  #   runPheno,
  #   seq(from = 1, to = nrow(runPheno))),
  #   recount::geo_characteristics)
  # # class(geochar.list)
  # # View(geochar.list)
  # # names(geochar.list)
  # # head(geochar.list)
  #
  #
  # ## collect column names
  # geo_colnames <- vector()
  # for (i in 1:length(geochar.list)) {
  #   geo_colnames <- union(geo_colnames, names((geochar.list[i])[[1]]))
  # }
  #
  # ## Convert the list to a data.frame
  # geochar.frame <- data.frame(matrix(nrow = length(geochar.list), ncol=length( geo_colnames))) ## Instantiate the data frame with the first sample
  # names(geochar.frame) <- geo_colnames
  # rownames(geochar.frame) <- names(geochar.list)
  # i <- 27
  # for (i in 1:length(geochar.list)) {
  #   sample.info <- (geochar.list[i])[[1]]
  #   geochar.frame[i, colnames(sample.info)] <- sample.info
  # }
  #
  # ## Replace column header "cells" by "cell.line"
  # ## This code comes from Leonardo, but we don't think we really need them for our purposes
  # ## The rbind does not work for some samples (e.g. SRP006574) because the number of columns of geochar varies from sample to sample.
  # # geochar <- do.call(rbind,
  # #                    lapply(geochar.list,
  # #                                  function(x) {
  # #                                    if ('cells' %in% colnames(x)) {
  # #                                      colnames(x)[colnames(x) == 'cells'] <- 'cell.line' ## Replace column header "cells" by "cell.line"
  # #                                      return(x)
  # #                                    } else {
  # #                                      return(x)
  # #                                    }
  # #                                  }))
  #
  # return(geochar.frame)
}
