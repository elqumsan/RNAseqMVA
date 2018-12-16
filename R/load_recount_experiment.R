#' @title  Loading one count table from Recount repository by given experiment ID and then having two option one made merge run to prevent
#' redundancy between multipl runs per sample.
#' @author Jacques van Helden and Mustafa Abuelqumsan
#' @description Load one count table from Recount for a given experiment ID, and optionally
#' merge the counts in order to avoid redundancy between multipl runs per sample.
#' @param recountID identifier of one study in ReCount database
#' @param parameters global and specific parameters for the analysis of this recountID
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
#' dim(recountData$result$runs$dataTable)
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
#' cor(log(recountData$result$runs$dataTable[, recountData$runPheno$geo_accession == "GSM1521620"]+1))
#'
#' ## Test correlation between 8 randomly selected biological replicates (distinct samples)
#' cor(log(recountData$merged$sampleCounts[, sample(1:ncol(recountData$merged$sampleCounts), size=8)]+1))
#'
#'
#' @return
#' A list containing: the count table, the pheno table, and some additional parameters (study ID, ...).
#'    \itemize{
#'         \item countsPerRun:  count table before any pre-proccessing procedure
#'         \item originalCounts: count table after merge runs proccess in order to get ride of multipl runs per sample.
#'    }
#'
#' @import recount
#' @import SummarizedExperiment
#' @import S4Vectors
#' @export
loadRecountExperiment <- function(recountID,
                                  parameters,
                                  # dir.workspace = parameters$dir$workspace,
                                  # mergeRuns = parameters$mergeRuns,
                                  # sampleIdColumn = parameters$sampleIdColumn, ## Alternative: use "sample"
                                  # classColumn = parameters$classColumn,
                                  # classColors = parameters$classColors,
                                  # variableType= parameters$variable.type,
                                  # verbose = parameters$verbose,
                                  forceDownload = FALSE,
                                  ...) {

#  LoadRequiredBioconductorPackages(c("recount", "SummarizedExperiment", "S4Vectors"))


  message.with.time("loadRecountExperiment()\trecountID = ", recountID)

  ## Check required parameters
  for (p in c("mergeRuns", "sampleIdColumn", "classColumn", "verbose", "studyPath")) {
    if (is.null(parameters[[p]])) {
      stop("Missing required parameter: '", p,
           "'.\n\tPlease check configuration file. ")
    } else {
      assign(p, parameters[[p]])
    }
  }

  ## Check required directory
  if (is.null(parameters$dir$workspace)) {
    stop("Missing required parameter: 'parameters$dir$workspace'.\n\tPlease check configuration file. ")
  } else {
    dir.workspace = parameters$dir$workspace
    dir.create(dir.workspace, recursive = TRUE, showWarnings = FALSE)
  }


  result <- list()

  #### Create studyPath directory ####
  if (!file.exists(parameters$studyPath)) {
    message("\tCreating directory to store Recount dataset ", recountID," in ", parameters$studyPath)
    dir.create(parameters$studyPath, recursive = TRUE, showWarnings = FALSE)
  }

  #### Define the file where the downloaded counts will be stored ####
  rseFile <- file.path(parameters$studyPath, "rse_gene.Rdata")
  parameters$rseFile <- rseFile

  #### Add parameters to the result ####
  result$parameters <- parameters

  #### Download the counts if required ####
  if ((forceDownload) || (!file.exists(rseFile))) {
    if (verbose) {
      message("\tDowloading counts from ReCount for study ", recountID)
    }
    url <- download_study(recountID, outdir = parameters$studyPath)
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
  rse <- scale_counts(rse_gene, by = "mapped_reads")

  #### Extract a matrix with the counts per feature for each run ####
  if (verbose) {
    message("\tExtracing table of counts per run")
  }


  ## Extract the count table
  dataTable <- assay(rse)
  if (verbose) {
    message("\tLoaded counts per run: ", nrow(dataTable), " features x ", ncol(dataTable), " runs.")
  }

  #### Extract pheno table ####
  ## The pheno table contains information about the columns of the RangedSeummaryExperiment.
  if (verbose) {
    message("\tBuilding pheno table")
  }
  phenoTable <- colData(rse) ## phenotype per run
#  geo.characteristics <- recount::geo_characteristics(phenoTable)
  geochar <- geocharFromPheno(runPheno = phenoTable)
  phenoTable <- cbind(phenoTable, geochar)
  # View(phenoTable)
  # names(phenoTable)


  countsPerRuns <- DataTableWithClasses(dataTable = dataTable,
                                         phenoTable = phenoTable,
                                         # classColumn = classColumn,
                                         # classColors = classColors,
                                         # variablesType = parameters$variables.type[1],
                                         dataType = "raw_counts_per_run",
                                         parameters = parameters)
  # class(countsPerRuns)
  summary(countsPerRuns)

  ################################################################
  ## Use either the merged or the original runs as original count table
  if (mergeRuns) {
    ## Merge runs if required
    if (verbose) { message("\tMerging run-wise counts by sample") }
    result$countsPerRun <- countsPerRuns
    result$originalCounts <- MergeRuns(runs = countsPerRuns) #,
                                  # classColumn  = classColumn ,
                                  # sampleIdColumn = sampleIdColumn,
                                  # verbose = verbose)
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
  # class(result$original) <- append(class(result$original), "DataTableWithClasses")

  message.with.time("Finished loading Recount experiment ID ", parameters$recountID)
  return(result)
}


## Mustafa; you make a method of StudyCase from thi rough code


## Histogram of mean counts per gene for a give class (Bone marrow here)
# counts <- as.vector(unlist(studyCases[[recountID]]$datasetsForTest$filtered$dataTable))
# hist(log2(counts + epsilon), breaks=100,
#      main = paste(recountID,
#                   " – Histogram of log2(counts)"),
#      xlab = "log2(counts)",
#      ylab = "Number of genes",
#      col="#CCBBFF")
# legend("topright",
#        legend = paste("Max counts per gene =",
#                       prettyNum(max(counts), big.mark = ",")))


# normcounts <-studyCases[[recountID]]$datasetsForTest$norm$dataTable
# gene.mean.per.class <- by(t(normcounts), INDICES = studyCases[[recountID]]$datasetsForTest$norm$classLabels, FUN = colMeans)
#
# epsilon <- 0.1
# x1 <- gene.mean.per.class$`Bone marrow` + epsilon
# x2 <- gene.mean.per.class$`Heparinised blood` + epsilon
#
#
# ## XY plot
# plot (x = log2(x1),
#       y = log2(x2),
#       main = paste(recountID, "\nlog2(scaled counts) per gene"),
#       xlab = "Bone marrow ",
#       ylab = "Heparinised blood",
#       las = 1,
#       col = densCols(x = log2(x1), y = log2(x2)),
#       panel.first = grid())
# abline(a = 0, b = 1, col = "black", lwd = 1)
#
# ## MA plot
# A <- (log2(x1) + log2(x2))/2
# M <- log2(x1) - log2(x2)
# plot (x = A, y = M,
#       main = paste(recountID, "\nMA plot"),
#       xlab = "A", ylab = "M",
#       col = densCols(x = A, y = M),
#       panel.first = grid())
# abline(h = 0, col = "black", lwd = 1)
#
#
# # plot(gene.mean.per.class$`Bone marrow` + epsilon,
# #      gene.mean.per.class$`Heparinised blood` + epsilon,
# #      log="xy")

