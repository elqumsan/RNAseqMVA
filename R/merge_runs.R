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

  unique.samples <- as.vector(unique(unlist(runPhenoTable[, sampleIdColumn])))
  message("\tMerging runs from count table (", nrow(countsPerRun), " features, ",
          ncol(countsPerRun), " runs), ",
          length(unique.samples), " unique samples. ")
  sample.nb <- length(unique.samples)
  sample.counts <- data.frame(matrix(nrow=nrow(countsPerRun),
                                     ncol=length(unique.samples)))
  names(sample.counts) <- unique.samples
  rownames(sample.counts) <- rownames(countsPerRun)
  # View(sample.counts)


  ## TO DO: see if we can use some apply or do.call() function
  ## to replace this loop, since "for" is an heresy in R
  s <- 0
  for (sample in unique.samples) {
    s <- s + 1
    # if (verbose) { message("\tmerging counts for sample ", s, "/", sample.nb, " ", sample) }
    runs <- grep(pattern = sample, x = runPhenoTable[, sampleIdColumn])
    if (length(runs ) > 1) {
      sample.counts[,sample] <- apply(countsPerRun[,runs],1,sum)
    } else {
      sample.counts[,sample] <- countsPerRun[,runs]
    }
  }

  ## Prepare a phenotable with all fields that are identical between runs
  ## This is tricky.
  sample.rows <- pmatch(unique.samples, unlist(runPhenoTable[sampleIdColumn]))
  sampleFields <- vector()
  for (field in names(runPhenoTable)) {
    nb.val <- length(unique(unlist(runPhenoTable[,field])))
    if (nb.val  <= sample.nb) {
      # if (verbose) { message("\tSample field ", field, "; values: ",  nb.val) }
      sampleFields <- append(sampleFields, field)
    }
  }
  samplePheno <- runPhenoTable[sample.rows, sampleFields]
  rownames(samplePheno) <- samplePheno[,sampleIdColumn]
  # View(samplePheno)
  sample.nb <- ncol(sample.counts)
  gene.nb <- nrow(sample.counts)

  message("\tCount table contains ", sample.nb, " samples (columns) and ", gene.nb, " genes (rows). ")
  result <- list(sampleCounts = sample.counts,
                 samplePheno = samplePheno,
                 sampleFields = sampleFields,
                 sampleIdColumn = sampleIdColumn)

  message.with.time("Finished MergeRuns()\t", "recountID = ", parameters$recountID)
  return(result)
}
