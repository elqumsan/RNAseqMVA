#' @title Merge read counts from different runs by sample.
#' @author Jacques van Helden
#' @description In ReCount data, some samples are sequenced over several runs, which increases the coverage but creates problems
#'  # because technical replicates come from the same sample and are thus not independent.
#' # We thus merge the runs per sample, in order to obtain a single vector of read counts per sample.
#' @param countsPerRun an object of the class CounttableWithPheno
#' @param sampleIdColumn="geo_accession"  name of the column of the pheno table which contains the sample IDs. GEO accession is preferred because it is widely used, but it is sometimes not defined. In such case, the "sample" column should be used.
#' @param verbose=FALSE if TRUE, write messages to indicate the progressing of the tasks
#' @export
MergeRuns <- function(runs,
                      sampleIdColumn = "geo_accession",
                      classColumn  = parameters$classColumn ,
                      verbose=FALSE) {
  message.with.time("MergeRuns()\t", "recountID = ", parameters$recountID)

  ## Check that the run argument belobgs to the right class
  class(runs)

  ## Extract gene names from the run count table
  geneNames <- rownames(runs$countTable)
  nbGenes <- length(runs$geneNames)

  ## Extract sample names from the specified column of the pheno table
  sampleNames <- as.vector(unique(unlist(runs$phenoTable[, sampleIdColumn])))
  nbSamples <- length(sampleNames)

  message("\tMerging runs from count table (",
          nbGenes, " genes, ",
          ncol(runs$countTable), " runs), ",
          nbSamples, " unique samples. ",
          "\n\tColumn from the pheno table used todefine sample IDs: ", sampleIdColumn)
  ## Build an empty count table with the right dimensions
  countTable <- data.frame(matrix(nrow=nrow(runs$countTable),
                                     ncol=nbSamples))
  colnames(countTable) <- sampleNames
  rownames(countTable) <- geneNames
  # View(countTable[1:10, 1:20])

  ## Collect the counts for each sample by summing the counts of the corresponding runs.
  s <- 0
  for (sample in sampleNames) {
    s <- s + 1
    # if (verbose) { message("\tmerging counts for sample ", s, "/", mergedRuns$nbSamples, " ", sample) }
    iterate <- grep(pattern = sample, x = runs$phenoTable[, sampleIdColumn])
    if (length( iterate ) > 1) {
      countTable[,sample] <- apply(runs$countTable[, iterate],1,sum)
    } else {
      countTable[,sample] <- runs$countTable[, iterate]
    }
  }
  # View(countTable[1:10, 1:20])

  ## Prepare a phenotable with all fields that are identical between runs
  ## This is tricky.
  sample.rows <- pmatch(sampleNames, unlist(runs$phenoTable[sampleIdColumn]))
  sampleFields <- vector()
  for (field in names(runs$phenoTable)) {
    nb.val <- length(unique(unlist(runs$phenoTable[,field])))
    if (nb.val  <= nbSamples) {
      # if (verbose) { message("\tSample field ", field, "; values: ",  nb.val) }
      sampleFields <- append(sampleFields, field)
    }
  }
  phenoTable <- runs$phenoTable[sample.rows, sampleFields]
  rownames(phenoTable) <- phenoTable[,sampleIdColumn]
  # View(mergedRuns$phenoTable)
  if (nbSamples!=  ncol(countTable)) {
    stop("The number of columns in the count table must equal the number of samples. ")
  }

#  message("\tCount table contains ", nbSamples, " samples (columns) and ", nbGenes, " genes (rows). ")

  ## Instantiate a new object of the class countTableWithClasses, with the merged runs
  mergedRuns <- countTableWithClasses(countTable = countTable,
                                    phenoTable = phenoTable,
                                    classColumn = classColumn,
                                    classesColors = runs$classesColors,
                                    variablesType = parameters$variables.type[1],
                                    dataType = "raw_counts_per_sample")
  # class(mergedRuns)
  if (verbose) {
    summary(mergedRuns)
  }

  message.with.time("Finished MergeRuns()\t", "recountID = ", parameters$recountID)


  return(mergedRuns)
}
