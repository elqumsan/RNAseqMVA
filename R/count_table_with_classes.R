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
#' @param ID=parameters$recountID identifier associated to the count table (by default, the RecountID, but can be specified with custom identifiers)
#' @param sampleNames=colnames(countTable) sample names (by default, automatically taken from the column names of the count table)
#' @param geneNames=townames(countTable) gene names (by default, automatically taken from the row names of the count table)
#' @param dataType="raw counts" data type, free text (e.g. raw counts, log2-transformed counts, log2 normalised counts  ...).
#'
#' @export

countTableWithClasses <- function(countTable,
                                  phenoTable,
                                  classColumn,
                                  ID = parameters$recountID,
                                  sampleNames = colnames(countTable),
                                  geneNames = rownames(countTable),
                                  dataType = "raw counts"
) {

  ## Built a list from the input parameters
message.with.time("\tCreate object has the countTableWithClasses attribute" )

  object <- structure(
    list(
      ID = ID,
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
#  message("\t\t\n giving the summary of the created object")
  cat("countTableWithClasses\n")
  cat("\tData type         \t", x$dataType, "\n")
  cat("\tNumber of genes   \t", x$nbGenes, "\n")
  cat("\tNumber of samples \t", x$nbSamples, "\n")
  cat("\tNumber of classes \t", x$nbClasses, "\n")
  cat("\tSamples per classes \n")
  print(x$samplesPerClass)
  cat("\n")
}

print.countTableWithClasses <- function(x) {
  summary.countTableWithClasses(x)
}

#' @title Export the different fields of an aobject in tab-separated values text files.
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @param self an object, which must belong to a compatible class
#' @export
exportTables <- function (self, ...) {
  message("Exporting object of class ", class(self), " to tables")
#  message("Looking for a function named ", paste("exportTables", class(self), sep="."))
  UseMethod("exportTables", self)

}

#' @title Export the different fields of an object of class countTableWithClasses in tab-separated values text files.
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @param self an object, which must belong to a compatible class
#' @param export.dir export directory
#' @param file.prefix file prefix to build the different tables
#' @param extension=".tsv" extension for the expodrted files (tab-separated value)
#' @param export.dir export directory
#' @param file.prefix file prefix to build the different tables
#' @param extension=".tsv" extension for the expodrted files (tab-separated value)
#' @export
exportTables.countTableWithClasses <- function (self,
                                                export.dir,
                                                file.prefix,
                                                extension=".tsv") {
  message("Exporting countTableWithClasses object ", self[["ID"]], " to tables")
  message("\tExport directory\t", export.dir)
  message("\tFile prefix\t", file.prefix)

  ############## Exporting the count table ####################
  count.file <- file.path(export.dir, paste(file.prefix, self[["ID"]], "_count_table", extension, sep = ""))
  message("\tExporting count table in TSV file\t", count.file)
  write.table(self$countTable, file = count.file, row.names = FALSE, quote=FALSE, sep = "\t")

  ############## Exporting the  Pheno Table ####################
  pheno.file <- file.path(export.dir, paste(file.prefix, self[["ID"]], "_pheno_table", extension, sep = ""))
  message("\tExporting  pheno table in TSV file\t", pheno.file)
  write.table(self$phenoTable, file = pheno.file, row.names = FALSE, quote=FALSE, sep = "\t")

  ############## Exporting the class labels ####################
  classLabel.file <- file.path(export.dir, paste(file.prefix, self[["ID"]], "_class_labels", extension, sep = ""))
  message("\tExporting class labels in TSV file\t", classLabel.file)
  write.table(data.frame(sampleName = self$sampleNames, classLabel = self$classLabels),
                         file = classLabel.file, row.names = FALSE, quote=FALSE, sep = "\t")

}
