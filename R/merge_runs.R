#' @title Merge read counts from different runs by sample.
#' @author Jacques van Helden
#' @description In ReCount data, some samples are sequenced over several runs, which increases the coverage but creates problems
#' # because technical replicates come from the same sample and are thus not independent.
#' # We thus merge the runs per sample, in order to obtain a single vector of read counts per sample.
#' @param countsPerRun a data frame containing read counts, one column per run and one row per gene
#' @param pheno.table a data frame containing the description of each run. This must include a column with the sample ID associated to each run
#' @param sampleIdColumn="geo_accession"  name of the column of the pheno table which contains the sample IDs. GEO accession is preferred because it is widely used, but it is sometimes not defined. In such case, the "sample" column should be used.
#' @param verbose=FALSE if TRUE, write messages to indicate the progressing of the tasks
#' @export
MergeRuns <- function(countsPerRun,
                      runPhenoTable,
                      sampleIdColumn = "geo_accession",
                      verbose=FALSE) {
  message.with.time("MergeRuns()\t", "recountID = ", parameters$recountID)
  result <- list()
  result$sampleIdColumn <- sampleIdColumn

  # dim(countsPerRun)

  result$sampleNames <- as.vector(unique(unlist(runPhenoTable[, sampleIdColumn])))
  result$nbSamples <- length(result$sampleNames)
  message("\tMerging runs from count table (",
          nrow(countsPerRun), " features, ",
          ncol(countsPerRun), " runs), ",
          result$nbSamples, " unique samples. ")
  result$countTable <- data.frame(matrix(nrow=nrow(countsPerRun),
                                     ncol=length(result$sampleNames)))
  names(result$countTable) <- result$sampleNames
  rownames(result$countTable) <- rownames(countsPerRun)
  # View(result$countTable)


  ## TO DO: see if we can use some apply or do.call() function
  ## to replace this loop, since "for" is an heresy in R
  s <- 0
  for (sample in result$sampleNames) {
    s <- s + 1
    # if (verbose) { message("\tmerging counts for sample ", s, "/", result$nbSamples, " ", sample) }
    runs <- grep(pattern = sample, x = runPhenoTable[, sampleIdColumn])
    if (length(runs ) > 1) {
      result$countTable[,sample] <- apply(countsPerRun[,runs],1,sum)
    } else {
      result$countTable[,sample] <- countsPerRun[,runs]
    }
  }
  # View(result$countTable)

  ## Prepare a phenotable with all fields that are identical between runs
  ## This is tricky.
  sample.rows <- pmatch(result$sampleNames, unlist(runPhenoTable[sampleIdColumn]))
  result$sampleFields <- vector()
  for (field in names(runPhenoTable)) {
    nb.val <- length(unique(unlist(runPhenoTable[,field])))
    if (nb.val  <= result$nbSamples) {
      # if (verbose) { message("\tSample field ", field, "; values: ",  nb.val) }
      result$sampleFields <- append(result$sampleFields, field)
    }
  }
  result$phenoTable <- runPhenoTable[sample.rows, result$sampleFields]
  rownames(result$phenoTable) <- result$phenoTable[,sampleIdColumn]
  # View(result$phenoTable)
  if (result$nbSamples!=  ncol(result$countTable)) {
    stop("The number of columns in the count table must equal the number of samples. ")
  }
  result$nbGenes <- nrow(result$countTable)

  message("\tCount table contains ", result$nbSamples, " samples (columns) and ", result$nbGenes, " genes (rows). ")

  message.with.time("Finished MergeRuns()\t", "recountID = ", parameters$recountID)

  ## S3 way of assigning a class of the result
  class(result) <- "CountTableWithDoc"
  class(result)
  names(result)

  return(result)
}
