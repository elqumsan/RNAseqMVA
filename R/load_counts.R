#' @title  Loading the count table, the Pheno table and vector of class labels
#' @author Mustafa ABUELQUMSAN and Jacques van Helden
#' @description Given a ReCount ID, load the count table, the pheno table
#' and the vector of class labels, prepares and filters the count table, via
#' wrapper runs by sample.
#' @param recountID  number of experiment
#' @param mergeRuns to wrapper the rus by sample
#' @param sampleIdColumn="geo_accession"  name of the column of the pheno table which contains the sample IDs.
#' @param classColumn=parameters$classColumn name of a column of the pheno table which contains the class labels.
#' In some cases classes must be built by concatenating several columns of the pheno table (e.g. "tissue" and "cell.type" for dataset SRP057196),
#' This can be achieved by providing a vector of column names from the pheno table. In this case, class names
## are built by concantenating the values in the specified columns (separated by "_").
#' @param minSamplesPerClass=parameters$minSamplesPerClass min number of samples per class for the filtering.  Passed to filterCountTable()
#' @param dir.workspace=parameters$dir$workspace Directory in which the data and results will be stored.
#' @param na.rm=TRUE if TRUE, omit all the rows (genes) that contain NA values. True by default because these NA values are problematic for classification methods.
#' @param ...  additional parameters are passed to loadRecountExperiment()
#'
#' @return
#'
#' @examples
#'
#' ##############################################
#' ## Load a RNA-seq dataset and merge counts per sample
#' x <- loadCounts( recountID = "SRP048759", mergeRuns=T, classColumn = "tissue")
#'
#' ################################################################
#' ## RNA-seq data combining 2 columns to define sample classes (classLabels vector)
#' ## Reduce the min number of samples per class for this dataset.
#' x <- loadCounts( recountID = "SRP057196", mergeRuns=T,
#'     classColumn = c("tissue", "cell.type"), minSamplesPerClass=5)
#'
#' @import recount
#' @import SummarizedExperiment
#' @import S4Vectors
#' @import caret
#' @export
loadCounts <- function(recountID = parameters$recountID,
                       classColumn = parameters$classColumn,
                       mergeRuns = parameters$mergeRuns,
                       sampleIdColumn = parameters$sampleIdColumn, ## Alternative: use "sample"
                       minSamplesPerClass = parameters$minSamplesPerClass,
                       dir.workspace = parameters$dir$workspace,
                       na.rm = parameters$na.rm,
                       ... ) {
  message.with.time("Starting loadCounts() for Recount ID ", recountID)

  ################################################
  # loading required libraries and install them if required
  # requiredpackage <- c("caret")
  # RequiredCRANPackages(requiredpackage)

  ################################################
  # loading count Data from recount_experiment, Via our wrapper which will Automatically merge the runs by
  # sample in order to obtain sample-wise count rather than run-wise counts.
  experiment <- loadRecountExperiment(recountID=recountID,
                                      mergeRuns=mergeRuns,
                                      sampleIdColumn=sampleIdColumn,
                                      dir.workspace = dir.workspace,
                                      classColumn = classColumn,
                                      na.rm = na.rm,
                                     ...)


  ## JhV: I suppress this na filtering from here because we want to keep the NA values in
  ## the original count table, and anyway wae do eliminate them in the filtering function.

  # ## If requested, suppress the NA values
  # if (na.rm) {
  #   message("\tSuppressing rows (genes) with NA values")
  #   countTable <- na.omit(experiment$originalCounts$countTable)
  #   # dim(countTable)
  # }

  ## This control is not necessary anymore since this is done when building the countTableWithClasses object.
  # ##### Check the dimensions of original experiment #####
  # if (experiment$originalCounts$nbSamples != length(experiment$originalCounts$classLabels)){
  #   stop("The number of samples (", experiment$originalCounts$nbSamples, ") differs from the number class labels (", length(experiment$original$classLabels),")")
  # }


  message("\toriginal countTable contains ",
          nrow(experiment$originalCounts$countTable), " rows (genes) and ",
          ncol(experiment$originalCounts$countTable), " columns (samples).")
  #countTable <-as.data.frame(t(experiment$runCounts))

  # ############## Exporting the original raw count data ####################
  # file.name <- file.path(tsv.dir, paste("original_counts_", parameters$recountID, ".tsv", sep = ""))
  # message("\tSaving original raw counts table in TSV file\t", file.name)
  # write.table(experiment$originalCounts$countTable, file = file.name, row.names = FALSE, quote=FALSE, sep = "\t")
  #
  # ############## Exporting the original Pheno Table ####################
  # file.name <- file.path(tsv.dir, paste("original_pheno_table_", parameters$recountID, ".tsv", sep = ""))
  # message("\tSaving original pheno table in TSV file\t", file.name)
  # write.table(experiment$originalCounts$phenoTable, file = file.name, row.names = FALSE, quote=FALSE, sep = "\t")
  #
  #
  # ##### Export a table with class labels + phenotypic information #####
  # message("\tSaving original pheno table + class labels information in TSV file\t", file.name)
  # characteristics.string <- unlist(lapply(experiment$originalCounts$phenoTable$characteristics, paste, collapse="; "))
  # pheno.data.frame <- data.frame(sample = experiment$originalCounts$phenoTable$sample,
  #                               # class = experiment$original$phenoTable,
  #                                description = characteristics.string)
  # # dim(pheno.data.frame)
  # # View(pheno.data.frame)
  #
  # file.name <- file.path(tsv.dir, paste("phenoDataframe_",parameters$recountID,".tsv", sep = ""))
  # message("\tExporting sample ID's and description in TSV file\t", file.name)
  # write.table(pheno.data.frame, file = file.name, quote=FALSE, row.names = FALSE, sep = "\t")


  ################################################
  #### Filter zero-variance and near-zero variance variables from the count table #####
  experiment$filtered <- filterCountTable(rawCounts = experiment$originalCounts,
                                          nearZeroVarFilter = parameters$nearZeroVarFilter,
                                          minSamplesPerClass = parameters$minSamplesPerClass)
  # class(experiment$filtered)
  # summary(experiment$filtered)
  # dim(countTable)


  ##### Check the dimensions of filtered experiment #####
  if (length(experiment$filtered$classLabels) != ncol(experiment$filtered$countTable)) {
    stop("invaled number of filtered class labels, (", length(experiment$filtered$classLabels) ,"). ",
         "must equal the number of individuals in Count Table (",ncol(experiment$filtered$countTable),").")
  }


  #### Build a list with the results of the loadCounts() function
  experiment$originalCounts$geo.characteristics <- geo_characteristics(experiment$originalCounts$phenoTable)
  experiment$filtered$geo.characteristics <- geo_characteristics(experiment$filtered$phenoTable)





  message.with.time("Finished loadCounts() for Recount experiment ID ", parameters$recountID)
  # experiment$countsPerRun <- NULL
  # experiment$runPhenoTable <- NULL
  # experiment$original <- NULL
  # experiment$filtered <- NULL
  # experiment$samples.per.class <- NULL
  # experiment$merged <- NULL
  # experiment$geo.characteristics <- NULL

  return(experiment)
}
