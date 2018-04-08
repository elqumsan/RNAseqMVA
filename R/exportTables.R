#' @title Export the different fields of an object in tab-separated values text files.
#' @description  Such target for ease the function for the biologists by converting the Rdata object for the tab-separated delimited file for facilitate the difficulties for extrating countTables and phenoTables
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @param self an object, which must belong to a compatible class
#' @export
exportTables <- function (self, ...) {
  message("Exporting object of class ", class(self), " to tables")
  #  message("Looking for a function named ", paste("exportTables", class(self), sep="."))
  UseMethod("exportTables", self)
}

#' @title default export table (for closing the NextMethod calls)
#' @author Jacques van Helden
#' @export
exportTables.default <- function (self, ...) {
  message("Finished exporting tables")
  #  message("Looking for a function named ", paste("exportTables", class(self), sep="."))
  return()
}



#' @title Export the different fields of an object of class countTableWithClasses in tab-separated values text files.
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @param self an object, which must belong to a compatible class
#' @param export.dir export directory
#' @param file.prefix file prefix to build the different tables
#' @param extension=".tsv" extension for the exported files (tab-separated values)
#' @export
exportTables.countTableWithClasses <- function (self,
                                                export.dir,
                                                file.prefix,
                                                extension=".tsv") {
  message("\tExporting countTableWithClasses object ", self$ID, " to tables")
  message("\t\tExport directory\t", export.dir)
  dir.create(export.dir, showWarnings = FALSE, recursive = TRUE)

  message("\t\tFile prefix\t", file.prefix)
  # UseMethod("exportTables", self)

  ############## Exporting the count table ####################
  count.file <- file.path(export.dir, paste(file.prefix, self$ID, "_count_table", extension, sep = ""))
  message("\t\tExporting count table in TSV file\t", count.file)
  write.table(self$countTable, file = count.file, row.names = FALSE, quote=FALSE, sep = "\t")

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

#' @title Export the tables specific to the class countTableWithTrainTestSets
#' @author Jacques van Helden
#' @param self an object, which must belong to class countTableWithTrainTestSets
#' @param export.dir export directory
#' @param file.prefix file prefix to build the different tables
#' @param extension=".tsv" extension for the exported files (tab-separated values)
#' @export
exportTables.countTableWithTrainTestSets <- function (self,
                                                export.dir,
                                                file.prefix,
                                                extension=".tsv") {
  message("\tExporting countTableWithTrainTestSets object ", self$ID, " to tables")
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




