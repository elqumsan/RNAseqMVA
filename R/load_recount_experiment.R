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
#' dim(recountData$result$runs$countTable)
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
#' cor(log(recountData$result$runs$countTable[, recountData$runPheno$geo_accession == "GSM1521620"]+1))
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
                                  classColumn = parameters$classColumn,
                                  classesColors = parameters$classesColors,
                                  variableType= parameters$variable.type,
                                  forceDownload = FALSE,
                                  verbose = parameters$verbose,
                                  ...) {
  message.with.time("loadRecountExperiment()\trecountID = ", recountID)

  result <- list()

  studyPath <- file.path(dir.workspace, "data", recountID)

  #### Create studyPath directory ####
  if (!file.exists(studyPath)) {
    message("\tCreating directory to store Recount dataset ", recountID," in ", studyPath)
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
      message("\tDowloading counts from ReCount for study ", recountID)
    }
    url <- download_study(recountID, outdir = studyPath)
    result$param["url"] <- url
  }


  #### Load in memory data from the recount database ####
  if (verbose) {
    message("\tLoading counts from local file ", rseFile)
  }
  load(rseFile)

  #### Scale counts by mapped reads, in order to to get read counts per gene ####
  if (verbose) {
    message("\tScaling counts")
  }
  rse <- scale_counts(rse_gene, by="mapped_reads")

  #### Extract a matrix with the counts per feature for each run ####
  if (verbose) {
    message("\tExtracing table of counts per run")
  }


  ## Extract the count table
  countTable <- assay(rse)
  if (verbose) {
    message("\tLoaded counts per run: ", nrow(countTable), " features x ", ncol(countTable), " runs.")
  }

  #### Extract pheno table ####
  ## The pheno table contains information about the columns of the RangedSeummaryExperiment.
  if (verbose) {
    message("\tBuilding pheno table")
  }
  phenoTable <- colData(rse) ## phenotype per run
  geo.characteristics <- recount::geo_characteristics(phenoTable)
  geochar <- geocharFromPheno(phenoTable)
  phenoTable <- cbind(phenoTable, geochar)
  # View(phenoTable)
  # names(phenoTable)


  countsPerRuns <- countTableWithClasses(countTable = countTable,
                                       phenoTable = phenoTable,
                                       classColumn = classColumn,
                                       classesColors = classesColors,
                                       variablesType = parameters$variables.type["all"],
                                       dataType = "raw_counts_per_run")
  # class(countsPerRuns)
  summary(countsPerRuns)

  ################################################################
  ## Use either the merged or the original runs as original count table
  if (mergeRuns) {
    ## Merge runs if required
    if (verbose) { message("\tMerging run-wise counts by sample") }
    result$countsPerRun <- countsPerRuns
    result$originalCounts <- MergeRuns(countsPerRuns,
                                  classColumn  = classColumn ,
                                  sampleIdColumn = sampleIdColumn,
                                  verbose = verbose)
  } else {
    result$originalCounts <- countsPerRun
  }
  # class(result$originalCounts)
  # summary(result$originalCounts)


  # ################################################################
  # ## Specify sample classes (classLabels) by extracting information about specified class columns
  # if (is.null(classColumn) || (length(classColumn) < 1)) {
  #   stop("classColumn must be defined. ")
  # } else if (length(classColumn) == 1) {
  #   result$original$classLabels <-  as.vector(result$original$phenoTable[, classColumn])
  # } else {
  #   ## Combine several columns to establish the classLabels
  #   result$original$classLabels <- apply(result$original$phenoTable[, classColumn], 1, paste, collapse="_")
  # }
  # table(classLabels)

  # result$original$classNames <- sort(unique(result$original$classLabels))
  # result$original$nbClasses <- length(result$original$classNames)
  # class(result$original) <- append(class(result$original), "CountTableWithClasses")

  message.with.time("Finished loading Recount experiment ID ", parameters$recountID)
  return(result)
}

