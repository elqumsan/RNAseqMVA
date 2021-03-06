#' @title Merge read counts from different runs by sample.
#' @author Jacques van Helden and Mustafa AbuELQumsan
#' @description In ReCount data, some samples are sequenced over several runs, which increases the coverage but creates problems
#' because technical replicates come from the same sample and are thus not independent.
#' We thus merge the runs per sample, in order to obtain a single vector of read counts per sample.
#' @param runs an object of the class DataTableWithClassess
#' @return  this is original cout table witout multipl runs per sample
#'  \itemize{
#'    \item mergedRuns that is the count table after marging the runs to be original count table after get ride of the multipl run per sample.
#'  }
#' @export
MergeRuns <- function(runs) {
                      # sampleIdColumn = "geo_accession",
                      # classColumn  = parameters$classColumn,
                      # verbose=FALSE) {
  message.with.time("MergeRuns()\t", "recountID = ", parameters$recountID)


  ## Check that the run argument belobgs to the class DataTableWithClassess
  # class(runs)
  ## Check the class of input object
  if (!is(runs, "DataTableWithClasses")) {
    stop("MergeRuns(): 'runs' parameter must belong to class DataTableWithClassess. ")
  }


  parameters <- runs$parameters
  for (p in c("sampleIdColumn", "classColumn", "verbose")) {
    if (is.null(parameters[[p]])) {
      stop("Missing required parameter: '", p,
           "'.\n\tPlease check configuration file. ")
    } else {
      assign(p, parameters[[p]])
    }
  }

  ## Extract gene names from the run count table
  featureNames <- rownames(runs$dataTable)
  nbGenes <- length(runs$featureNames)

  ## Extract sample names from the specified column of the pheno table
  sampleNames <- as.vector(unique(unlist(runs$phenoTable[, sampleIdColumn])))
  nbSamples <- length(sampleNames)

  message("\tMerging runs from count table (",
          nbGenes, " genes, ",
          ncol(runs$dataTable), " runs), ",
          nbSamples, " unique samples. ",
          "\n\tColumn from the pheno table used to define sample IDs: ", sampleIdColumn)
  ## Build an empty count table with the right dimensions
  dataTable <- data.frame(matrix(nrow=nrow(runs$dataTable),
                                     ncol=nbSamples))
  colnames(dataTable) <- sampleNames
  rownames(dataTable) <- featureNames
  # View(dataTable[1:10, 1:20])

  ## Collect the counts for each sample by summing the counts of the corresponding runs.
  s <- 0
  for (sample in sampleNames) {
    s <- s + 1
    # if (verbose) { message("\tmerging counts for sample ", s, "/", mergedRuns$nbSamples, " ", sample) }
    iterate <- grep(pattern = sample, x = runs$phenoTable[, sampleIdColumn])
    if (length( iterate ) > 1) {
      dataTable[,sample] <- apply(runs$dataTable[, iterate],1,sum)
    } else {
      dataTable[,sample] <- runs$dataTable[, iterate]
    }
  }
  # View(dataTable[1:10, 1:20])

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
  if (nbSamples !=  ncol(dataTable)) {
    stop("The number of columns in the count table must equal the number of samples. ")
  }

#  message("\tCount table contains ", nbSamples, " samples (columns) and ", nbGenes, " genes (rows). ")

  ## Instantiate a new object of the class DataTableWithClassess, with the merged runs
  mergedRuns <- DataTableWithClasses(dataTable = dataTable,
                                      phenoTable = phenoTable,
                                      # classColumn = classColumn,
                                      # classColors = runs$classColors,
                                      # variablesType = parameters$variables.type[1],
                                      dataType = "raw_counts_per_sample",
                                      parameters = runs$parameters)
  # class(mergedRuns)
  if (verbose) {
    summary(mergedRuns)
  }

  message.with.time("Finished MergeRuns()\t", "recountID = ", parameters$recountID)


  return(mergedRuns)
}
