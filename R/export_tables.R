#' @title Export the different fields of an object in tab-separated values text files.
#' @description  Such target for ease the function for the biologists by converting the Rdata object for the tab-separated delimited file for facilitate the difficulties for extrating dataTables and phenoTables
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @param self an object, which must belong to a compatible class
#' @export
exportTables <- function (self, ...) {
  message("Exporting object of class ", class(self), " to tables")
  #  message("Looking for a function named ", paste("exportTables", class(self), sep="."))
  UseMethod("exportTables", self)
}

#' @title default export table (for closing the NextMethod calls)
#' @author Jacques van Helden and Mustafa AbuELQumsan.
#' @description  to print message which is inform us the exporting processes has been done.
#' @param self  an object, which must belong to a compatible class
#'
#' @export
exportTables.default <- function (self, ...) {
  message("Finished exporting tables")
  #  message("Looking for a function named ", paste("exportTables", class(self), sep="."))
  cat("\n") ## This avoids to print NULL in the console
  return()
}



#' @title Export the different fields of an object of class DataTableWithClasses in tab-separated values text files.
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @param self an object, which must belong to a compatible class
#' @param export.dir export directory
#' @param file.prefix file prefix to build the different tables
#' @param extension=".tsv" extension for the exported files (tab-separated values)
#' @export
exportTables.DataTableWithClasses <- function (self,
                                                export.dir = self$parameters$dir$tsv,
                                                file.prefix,
                                                extension=".tsv") {
  message("\tExporting DataTableWithClasses object ", self$ID, " to tables")
  message("\t\tExport directory\t", export.dir)
  dir.create(export.dir, showWarnings = FALSE, recursive = TRUE)

  message("\t\tFile prefix\t", file.prefix)
  # UseMethod("exportTables", self)

  ############## Exporting the count table ####################
  count.file <- file.path(export.dir, paste(file.prefix, self$ID, "_count_table", extension, sep = ""))
  message("\t\tExporting count table in TSV file\t", count.file)
  write.table(self$dataTable, file = count.file, row.names = FALSE, quote=FALSE, sep = "\t")

  ############## Exporting the  Pheno Table ####################
  pheno.file <- file.path(export.dir, paste(file.prefix, self$ID, "_pheno_table", extension, sep = ""))
  message("\t\tExporting  pheno table in TSV file\t", pheno.file)
  write.table(self$phenoTable, file = pheno.file, row.names = FALSE, quote=FALSE, sep = "\t")

  ############## Exporting the class labels ####################
  classLabel.file <- file.path(export.dir, paste(file.prefix, self$ID, "_class_labels", extension, sep = ""))
  message("\t\tExporting class labels in TSV file\t", classLabel.file)
  write.table(data.frame(sampleName = self$sampleNames, classLabel = self$classLabels),
              file = classLabel.file, row.names = FALSE, quote=FALSE, sep = "\t")

  NextMethod("exportTables", self)
}

#' @title Export the tables specific to the class DataTableWithTrainTestSets
#' @author Mustafa AbuElQumsan and Jacques van Helden.
#' @param self an object, which must belong to class DataTableWithTrainTestSets
#' @param export.dir export directory
#' @param file.prefix file prefix to build the different tables
#' @param extension=".tsv" extension for the exported files (tab-separated values)
#' @export
exportTables.DataTableWithTrainTestSets <- function (self,
                                                export.dir,
                                                file.prefix,
                                                extension=".tsv") {
  message("\tExporting DataTableWithTrainTestSets object ", self$ID, " to tables")
  message("\t\tExport directory\t", export.dir)
  dir.create(export.dir, showWarnings = FALSE, recursive = TRUE)

  message("\t\tFile prefix\t", file.prefix)
  # UseMethod("exportTables", self)

  ############## Exporting train and test indices ####################
  trainIndices.file <- file.path(export.dir, paste(file.prefix, self$ID, "_trainindices", extension, sep = ""))
  message("\t\tExporting count table in TSV file\t", trainIndices.file)
  write.table(as.data.frame(self$trainTestProperties$trainIndices), file = trainIndices.file, row.names = FALSE, col.names = FALSE, quote=FALSE, sep = "\t")

  testIndices.file <- file.path(export.dir, paste(file.prefix, self$ID, "_testindices", extension, sep = ""))
  message("\t\tExporting count table in TSV file\t", testIndices.file)
  write.table(as.data.frame(self$testTestProperties$testIndices), file = testIndices.file, row.names = FALSE, col.names = FALSE, quote=FALSE, sep = "\t")


  NextMethod("exportTables", self)
}


#' @title Export the tables specific to the class StudyCase
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @description this calls exportTable on each of the datasets  (which belong to DataTableWithClasses and DataTableWithTrainTestSets)
#' @param self an object, which must belong to class StudyCase
#' @param export.dir export directory
#' @param file.prefix file prefix to build the different tables
#' @param extension=".tsv" extension for the exported files (tab-separated values)
#' @export
exportTables.StudyCase <- function(self,
                                  export.dir = self$parameters$dir$tsv,
                                  extension=".tsv") {

  message("\tExporting StudyCase object ", self$ID, " to tables")
  message("\t\tExport directory\t", export.dir)

  ## Export raw counts per run
  exportTables(self$rawData$countsPerRun,
               export.dir = file.path(export.dir, self$ID),
               file.prefix = "counts_per_run_")

  ## Export raw counts per sample
  exportTables(self$rawData$countsPerSample,
               export.dir = file.path(export.dir, self$ID),
               file.prefix = "counts_per_sample_")

  ## Export filtere counts
  exportTables(self$datasetsForTest$filtered,
               export.dir = file.path(export.dir, self$ID),
               file.prefix = "filtered_counts_")

  ## Export normalized counts
  exportTables(self$datasetsForTest$norm,
               export.dir = file.path(export.dir, self$ID),
               file.prefix = "norm_counts_")

  ## Export log2-transformed normalised counts
  exportTables(self$datasetsForTest$log2norm,
               export.dir = paste(export.dir, self$ID, sep = "/"),
               file.prefix = "log2norm_counts_")

  ## Export principal components of log2-transformed normalised counts
  exportTables(self$datasetsForTest$log2normPCs,
               export.dir = file.path(export.dir, self$ID),
               file.prefix = "log2norm_counts_")

  NextMethod("exportTables", self)

}

