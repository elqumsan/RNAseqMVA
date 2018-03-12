#' @title constructor of the countTableWithClasses class
#' @author Mustafa AbuElQumsan and Jacques van helden
#' @description This class contains count tables
#' (typically used to represent expression profiles in RNA-seq analysis)
#' documented with a Pheno table (description of each sample) and completed
#' with some additional fields for convenience: sample names, gene names, ...
#' @param countTable a data.frame with one row per feature (e.g. gene) and
#' one column per sample.
#' @param phenoTable a data frame describing each sample: one row per sample and one column per attribute of a sample.
#' @param classColumn a vector indicating one or several columns of the pheno table that will define the sample class labels. If umtiple columns are specified, they will be concatenated to build sample labels
#' @param sampleNames=colnames(countTable) sample names (by default, automatically taken from the column names of the count table)
#' @param geneNames=townames(countTable) gene names (by default, automatically taken from the row names of the count table)
#' @param dataType="raw counts" data type, free text (e.g. raw counts, log2-transformed counts, log2 normalised counts  ...).
#'
#' @export

countTableWithClasses <- function(countTable,
                                phenoTable,
                                classColumn,
                                sampleNames = colnames(countTable),
                                geneNames = rownames(countTable),
                                dataType = "raw counts"
) {

  ## Built a list from the input parameters
message.with.time("\tCreate object has the countTableWithClasses attribute" )

  object <- structure(
    list(
      countTable = countTable,
      phenoTable = phenoTable,
      sampleNames = sampleNames,
      nbSamples = ncol(countTable),
      geneNames = geneNames,
      nbGenes = nrow(countTable),
      dataType = dataType
    ),
    class="countTableWithClasses")
  # names(object)
  # class(object)
  # attributes(object)


  #### Check consistency of input parameters ####

  ## Check rows of pheno table
  if (nrow(object$phenoTable) != ncol(object$countTable)) {
    stop("countTableWithClasses(): inconsistent dimensions of phenoTable (",
         nrow(phenoTable), " rows) and countTable (",
         ncol(countTable), " columns).")
  }

  ## Check sample names
  if (length(object$sampleNames) != ncol(countTable)) {
    stop("countTableWithClasses(): inconsistent dimensions of sampleNames (",
         length(object$sampleNames), " names) and countTable (",
         ncol(countTable), " columns).")
  }

  ## Check gene names
  if (length(object$geneNames) != nrow(countTable)) {
    stop("countTableWithClasses(): inconsistent dimensions of geneNames (",
         length(object$geneNames), " names) and countTable (",
         nrow(countTable), " columns).")
  }


  ################################################################
  ## Specify sample classes (classLabels) by extracting information about specified class columns
  if (is.null(classColumn) || (length(classColumn) < 1)) {
    stop("classColumn must be defined. ")
  } else if (length(classColumn) == 1) {
    object$classLabels <-  as.vector(object$phenoTable[, classColumn])
  } else {
    ## Combine several columns to establish the classLabels
    object$classLabels <- apply(object$phenoTable[, classColumn], 1, paste, collapse="_")
  }
  # table(classLabels)
  object$classNames <- unique(sort(object$classLabels))
  object$nbClasses <- length(object$classNames)
  object$samplesPerClass <- (as.data.frame.table(table(object$classLabels)))
  colnames(object$samplesPerClass) <- c("Class", "nbSamples")

  message("\t\tfinishing from creating the object with countTablewithClasses attribute")
  return(object)
}


summary.countTableWithClasses <- function(x) {
  cat("countTableWithClasses\n")
  cat("\tData type         \t", x$dataType, "\n")
  cat("\tNumber of genes   \t", x$nbGenes, "\n")
  cat("\tNumber of samples \t", x$nbSamples, "\n")
  cat("\tNumber of classes \t", x$nbClasses, "\n")
  cat("\tSamples per classes \n")
  print(x$samplesPerClass)
  cat("\n")
}

message("\t\t\n giving the summary of the created object")
print.countTableWithClasses <- function(x) {
  summary.countTableWithClasses(x)
}
