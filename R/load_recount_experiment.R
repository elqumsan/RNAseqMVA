#' @title Extract the conditions from the "characteristics" column of the coldata.
#' @description This is a bit tricky: we have to parse a string describing several attributes.
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


#' @title Load one count table from ReCount and optionally merge the counts per sample.
#' @author: Jacques van Helden and Mustafa Abuelqumsan
#' @description Load one count table from Recount for a given experiment ID, and optionally
#' merge the counts in order to avoid redundancy between multipl runs per sample.
#' @param recountID=parameters$recountID identifier of one study in ReCount database
#' @param dir.workspace=parameters$dir$workspace path of the folder to store the results (in one separate sub-directory per recountID)
#' @param mergeRuns=parameters$mergeRuns if TRUE, read counts will be merged for each sample
#' @param sampleIdColumn=parameters$sampleIdColumn  name of the column of the pheno table which contains the sample IDs.
#' This information is passed to MergeRuns().
#' @param verbose=TRUE if TRUE, write messages to indicate the progressing of the tasks
#' @param forceDownload=FALSE by default, the data is downloaded only if it is not found in the studyPath folder.
#' If forceDownload is TRUE, the data will be downloaded irrespective of existing files.
#'
#' @examples
#'
#' ## Load counts for a Leukemia dataset
#' recountID <- "SRP048759"
#' recountData <- loadRecountExperiment(recountID)
#'
#' ## Check the dimension of the table with counts per run
#' dim(recountData$countsPerRun)
#'
#' ## Check the number of runs per sample
#' table(recountData$runPheno$geo_accession)
#'
#' ## Count the number of unique sample IDs in the run-wise description table (runPheno)
#' message("Number of runs = ", length(recountData$runPheno$geo_accession))
#' message("Number unique geo_accession values: ", length(unique(sort(recountData$runPheno$geo_accession))))
#'
#' ## Check the dimension of the table with counts per sample
#' dim(recountData$merged$sampleCounts)
#'
#' cor(log(recountData$countsPerRun[, recountData$runPheno$geo_accession == "GSM1521620"]+1))
#'
#' ## Test correlation between 8 randomly selected biological replicates (distinct samples)
#' cor(log(recountData$merged$sampleCounts[, sample(1:ncol(recountData$merged$sampleCounts), size=8)]+1))
#'
#'
#' @return
#' A list containing: the count table, the pheno table, and some additional parameters (study ID, ...).
#'
#'
#' @import recount
#' @import SummarizedExperiment
#' @import S4Vectors
#'
#' @export
loadRecountExperiment <- function(recountID = parameters$recountID,
                                  dir.workspace = parameters$dir$workspace,
                                  mergeRuns = parameters$mergeRuns,
                                  sampleIdColumn = parameters$sampleIdColumn, ## Alternative: use "sample"
                                  forceDownload = FALSE,
                                  verbose = parameters$verbose,
                                  ...) {
  result <- list()
  studyPath <- file.path(dir.workspace, "data", recountID)

  #### Create studyPath directory ####
  if (!file.exists(studyPath)) {
    message("Creating directory to store Recount dataset ", recountID," in ", studyPath)
    dir.create(studyPath, recursive = TRUE, showWarnings = FALSE)
  }

  #### Define the file where the downloaded counts will be stored ####
  rseFile <- file.path(studyPath, "rse_gene.Rdata")

  #### Add parameters to the result ####
  result$param <-
    c("recountID" = recountID,
      "studyPath" = studyPath,
      "mergeRuns" = mergeRuns,
      "rseFile" = rseFile)

  #### Download the counts if required ####
  if ((forceDownload == TRUE) || (!file.exists(rseFile))) {
    if (verbose) {
      message("Dowloading counts from ReCount for study ", recountID)
    }
    url <- download_study(recountID, outdir = studyPath)
    result$param["url"] <- url
  }


  #### Load in memory data from the recount database ####
  if (verbose) {
    message("Loading counts from local file ", rseFile)
  }
  load(rseFile)

  #### Scale counts by mapped reads, in order to to get read counts per gene ####
  if (verbose) {
    message("Scaling counts")
  }
  rse <- scale_counts(rse_gene, by="mapped_reads")
  result$rse <- rse

  #### Extract a matrix with the counts per feature for each run ####
  if (verbose) {
    message("Extracing table of counts per run")
  }
  countsPerRun <- assay(rse)
  # View(countsPerRun)
  result$countsPerRun <- countsPerRun
  if (verbose) {
    message("Loaded counts per run: ", nrow(countsPerRun), " features x ", ncol(countsPerRun), " runs.")
  }

  #### Extract pheno table ####
  ## The pheno table contains information about the columns of the RangedSeummaryExperiment.
  if (verbose) {
    message("Building pheno table")
  }
  runPheno <- colData(rse) ## phenotype per run
  # class(runPheno)
  # pheno <- runPheno ## A TRICK
  # View(runPheno)
  # names(runPheno)

  geochar <- geocharFromPheno(runPheno)
  runPhenoTable <- cbind(runPheno, geochar)

  # View(as.data.frame(geochar))

  # View(geochar)
  # head(geochar)

  # Build a pheno table with selected columns from coldata + the geodata we just extracted
  # runPhenoTable <- cbind(
  #   runPheno[, grep(pattern="(characteristics|sharq)", x=names(runPheno), invert=TRUE)],
  #   geochar)
  # View(phenoTable)
  # class(phenoTable)

  # ## MUSTAFA STOPE THE TWO FOLLOWING PARTS DUE TO THE error duplicate in row.mane

  # ## Extract a phenoTable with selected fields from the runPheno object
  # runPhenoTable2 <- data.frame(
  #   project = runPheno$project,
  #   sample = runPheno$sample,
  #   experiment = runPheno$experiment,
  #   run = runPheno$run,
  #   geo_accession = runPheno$geo_accession,
  #   characteristics = runPheno$characteristics@unlistData
  # )
  #
  # ## Missing: parse sub-fields from the "characteristics" field (somehow tricky)
  #
  # rownames(runPhenoTable2) <- runPhenoTable2$run
  # # View(runPhenoTable2)
  # # class(runPhenoTable2)
  # result$runPhenoTable2 <- runPhenoTable2

  # STOP Mustafa
  if (mergeRuns) {
    if (verbose) { message("Merging run-wise counts by sample") }
    result$merged <- MergeRuns(countsPerRun,
                               runPhenoTable,
                               sampleIdColumn = sampleIdColumn,
                               verbose = verbose)
  }

  message.with.time("Finished loading Recount experiment ID ", parameters$recountID)
  return(result)
}

