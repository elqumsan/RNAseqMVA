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

  ## S3 way of assigning a class of the mergedRuns
  mergedRuns <- list()
  class(mergedRuns) <- "CountTableWithPheno"
  mergedRuns$dataType <- "counts per sample"
  mergedRuns$sampleIdColumn <- sampleIdColumn

  # dim(countsPerRun)

  mergedRuns$geneNames <- rownames(countsPerRun)
  mergedRuns$nbGenes <- length(mergedRuns$geneNames)

  mergedRuns$sampleNames <- as.vector(unique(unlist(runPhenoTable[, sampleIdColumn])))
  mergedRuns$nbSamples <- length(mergedRuns$sampleNames)

  message("\tMerging runs from count table (",
          mergedRuns$nbGenes, " genes, ",
          ncol(countsPerRun), " runs), ",
          mergedRuns$nbSamples, " unique samples. ")
  mergedRuns$countTable <- data.frame(matrix(nrow=nrow(countsPerRun),
                                     ncol=length(mergedRuns$sampleNames)))
  names(mergedRuns$countTable) <- mergedRuns$sampleNames
  rownames(mergedRuns$countTable) <- rownames(countsPerRun)
  # View(mergedRuns$countTable)


  ## TO DO: see if we can use some apply or do.call() function
  ## to replace this loop, since "for" is an heresy in R
  s <- 0
  for (sample in mergedRuns$sampleNames) {
    s <- s + 1
    # if (verbose) { message("\tmerging counts for sample ", s, "/", mergedRuns$nbSamples, " ", sample) }
    runs <- grep(pattern = sample, x = runPhenoTable[, sampleIdColumn])
    if (length(runs ) > 1) {
      mergedRuns$countTable[,sample] <- apply(countsPerRun[,runs],1,sum)
    } else {
      mergedRuns$countTable[,sample] <- countsPerRun[,runs]
    }
  }
  # View(mergedRuns$countTable)

  ## Prepare a phenotable with all fields that are identical between runs
  ## This is tricky.
  sample.rows <- pmatch(mergedRuns$sampleNames, unlist(runPhenoTable[sampleIdColumn]))
  mergedRuns$sampleFields <- vector()
  for (field in names(runPhenoTable)) {
    nb.val <- length(unique(unlist(runPhenoTable[,field])))
    if (nb.val  <= mergedRuns$nbSamples) {
      # if (verbose) { message("\tSample field ", field, "; values: ",  nb.val) }
      mergedRuns$sampleFields <- append(mergedRuns$sampleFields, field)
    }
  }
  mergedRuns$phenoTable <- runPhenoTable[sample.rows, mergedRuns$sampleFields]
  rownames(mergedRuns$phenoTable) <- mergedRuns$phenoTable[,sampleIdColumn]
  # View(mergedRuns$phenoTable)
  if (mergedRuns$nbSamples!=  ncol(mergedRuns$countTable)) {
    stop("The number of columns in the count table must equal the number of samples. ")
  }

  message("\tCount table contains ", mergedRuns$nbSamples, " samples (columns) and ", mergedRuns$nbGenes, " genes (rows). ")

  message.with.time("Finished MergeRuns()\t", "recountID = ", parameters$recountID)


  return(mergedRuns)
}
