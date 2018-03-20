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




#' @title Export the different fields of an object of class countTableWithClasses in tab-separated values text files.
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @param self an object, which must belong to a compatible class
#' @param export.dir export directory
#' @param file.prefix file prefix to build the different tables
#' @param extension=".tsv" extension for the expodrted files (tab-separated value)
#' @export
exportTables.countTableWithClasses <- function (self,
                                                export.dir,
                                                file.prefix,
                                                extension=".tsv") {
  message("\tExporting countTableWithClasses object ", self[["ID"]], " to tables")
  message("\t\tExport directory\t", export.dir)
  message("\t\tFile prefix\t", file.prefix)

  ############## Exporting the count table ####################
  count.file <- file.path(export.dir, paste(file.prefix, self[["ID"]], "_count_table", extension, sep = ""))
  message("\t\tExporting count table in TSV file\t", count.file)
  write.table(self$countTable, file = count.file, row.names = FALSE, quote=FALSE, sep = "\t")

  ############## Exporting the  Pheno Table ####################
  pheno.file <- file.path(export.dir, paste(file.prefix, self[["ID"]], "_pheno_table", extension, sep = ""))
  message("\t\tExporting  pheno table in TSV file\t", pheno.file)
  write.table(self$phenoTable, file = pheno.file, row.names = FALSE, quote=FALSE, sep = "\t")

  ############## Exporting the class labels ####################
  classLabel.file <- file.path(export.dir, paste(file.prefix, self[["ID"]], "_class_labels", extension, sep = ""))
  message("\t\tExporting class labels in TSV file\t", classLabel.file)
  write.table(data.frame(sampleName = self$sampleNames, classLabel = self$classLabels),
              file = classLabel.file, row.names = FALSE, quote=FALSE, sep = "\t")

}

