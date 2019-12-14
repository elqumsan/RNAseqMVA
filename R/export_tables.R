#' @title Export the different fields of an object in tab-separated values text files.
#' @description  Such target for ease the function for the biologists by converting the Rdata object for the tab-separated delimited file for facilitate the difficulties for extrating dataTables and phenoTables
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @param self an object, which must belong to a compatible class
#' @export
exportTables <- function(self, ...) {
  message("\t\tExporting object of class ", class(self), " to tables")
  #  message("Looking for a function named ", paste("exportTables", class(self), sep="."))
  UseMethod("exportTables", self)
}

#' @title default export table (for closing the NextMethod calls)
#' @author Jacques van Helden and Mustafa AbuELQumsan.
#' @description  to print message which is inform us the exporting processes has been done.
#' @param self  an object, which must belong to a compatible class
#'
#' @export
exportTables.default <- function(self, ...) {
  message("\t\tFinished exporting tables")
  #  message("Looking for a function named ", paste("exportTables", class(self), sep="."))
  cat("\n") ## This avoids to print NULL in the console
 # return()
}



#' @title Export the different fields of an object of class DataTableWithClasses in tab-separated values text files.
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @param self an object, which must belong to a compatible class
#' @param export.dir export directory
#' @param file.prefix file prefix to build the different tables
#' @param extension=".tsv" extension for the exported files (tab-separated values)
#' @export
exportTables.DataTableWithClasses <- function(self,
                                                export.dir = self$parameters$dir$tsv,
                                                file.prefix,
                                                extension=".tsv") {
  message("\tExporting DataTableWithClasses object ", self$ID, " to tables")
  message("\t\tExport directory\t", export.dir)
  dir.create(export.dir, showWarnings = FALSE, recursive = TRUE)

  message("\t\tFile prefix\t", file.prefix)
  # UseMethod("exportTables", self)

  ############## Exporting the count table ####################
  count.file <- file.path(export.dir,
                          paste0(file.prefix, self$ID, "_", self$parameters$feature, "_count_table", extension))
  message("\t\tExporting count table to TSV file\t", count.file)
  write.table(self$dataTable, file = count.file, row.names = FALSE, quote = FALSE, sep = "\t")

  ############## Exporting the  Pheno Table ####################
  pheno.file <- file.path(export.dir, paste0(file.prefix, self$ID, "_", self$parameters$feature, "_pheno_table", extension))
  message("\t\tExporting  pheno table to TSV file\t", pheno.file)
  write.table(self$phenoTable, file = pheno.file, row.names = FALSE, quote = FALSE, sep = "\t")

  ############## Exporting the class labels ####################
  classLabel.file <- file.path(export.dir, paste0(file.prefix, self$ID, "_", self$parameters$feature, "_class_labels", extension))
  message("\t\tExporting class labels to TSV file\t", classLabel.file)
  write.table(data.frame(sampleName = self$sampleNames, classLabel = self$classLabels),
              file = classLabel.file, row.names = FALSE, quote = FALSE, sep = "\t")

  NextMethod("exportTables", self)
}

#' @title Export the tables specific to the class DataTableWithTrainTestSets
#' @author Mustafa AbuElQumsan and Jacques van Helden.
#' @param self an object, which must belong to class DataTableWithTrainTestSets
#' @param export.dir export directory
#' @param file.prefix file prefix to build the different tables
#' @param extension=".tsv" extension for the exported files (tab-separated values)
#' @export
exportTables.DataTableWithTrainTestSets <- function(self,
                                                export.dir,
                                                file.prefix,
                                                extension=".tsv") {
  message("\tExporting DataTableWithTrainTestSets object ", self$ID, " to tables")
  message("\t\tExport directory\t", export.dir)
  dir.create(export.dir, showWarnings = FALSE, recursive = TRUE)

  message("\t\tFile prefix\t", file.prefix)
  # UseMethod("exportTables", self)

  ############## Exporting train and test indices ####################
  trainIndices.file <- file.path(export.dir, paste0(file.prefix, self$ID, "_", self$parameters$feature,  "_trainindices", extension))
  message("\t\tExporting count table to TSV file\t", trainIndices.file)
  write.table(as.data.frame(self$trainTestProperties$trainIndices), file = trainIndices.file, row.names = TRUE, col.names = NA,  quote = FALSE, sep = "\t")

  testIndices.file <- file.path(export.dir, paste0(file.prefix, self$ID, "_", self$parameters$feature,"_testindices", extension))
  message("\t\tExporting count table to TSV file\t", testIndices.file)
  write.table(as.data.frame(self$testTestProperties$testIndices), file = testIndices.file, row.names = TRUE, col.names = NA,  quote = FALSE, sep = "\t")


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
                                  # file.prefix,
                                  extension=".tsv") {

  message("\tExporting StudyCase object ", self$ID, " to tables")
  message("\t\tExport directory\t", export.dir)

  ## Export raw counts per run
  ############## Exporting train and test indices ####################
  # exportTables(self$rawData$countsPerRun,
  #              export.dir = file.path(export.dir, self$ID),
  #              file.prefix = "counts_per_run_")

  countsPerRun.file <- file.path(export.dir, paste0(self$ID, "_", self$parameters$feature, "_counts_per_run", extension))
  message("\t\tExporting counts per run to TSV file\t", countsPerRun.file)
  write.table(as.data.frame(self$rawData$countsPerRun$dataTable), file = countsPerRun.file, row.names = TRUE, col.names = NA, quote = FALSE, sep = "\t")
  # system(paste("head -n 10 ", countsPerRun.file, " | cut -f 1-5"))
  # View(head(self$rawData$countsPerRun$dataTable))


  ## Export raw counts per sample
  countsPerSample.file <- file.path(export.dir, paste0(self$ID, "_", self$parameters$feature, "_counts_per_sample", extension))
  message("\t\tExporting count per sample to TSV file\t", countsPerSample.file)
  write.table(as.data.frame(self$rawData$countsPerSample$dataTable), file = countsPerSample.file, row.names = TRUE, col.names = NA,  quote = FALSE, sep = "\t")

  ## Export filtere counts
  # exportTables(self$datasetsForTest$filtered,
  #              export.dir = file.path(export.dir, self$ID),
  #              file.prefix = "filtered_counts_")
  filtered.file <- file.path(export.dir, paste0(self$ID, "_", self$parameters$feature, "_filtered", extension))
  message("\t\tExporting filtered DataTable to TSV file\t", filtered.file)
  write.table(as.data.frame(self$datasetsForTest$filtered$dataTable), file = filtered.file, row.names = TRUE, col.names = NA,  quote = FALSE, sep = "\t")

  ## Export the different tables of normalized values
  for (datasetName in names(self$datasetsForTest)) {
    file <- file.path(export.dir, paste0(self$ID, "_", self$parameters$feature, "_", datasetName, extension))
    message("\t\tExporting scaled DataTable to TSV file\t", file)
    write.table(as.data.frame(self$datasetsForTest[[datasetName]]$dataTable),
                file = file, row.names = TRUE, col.names = NA,  quote = FALSE, sep = "\t")

  }

  NextMethod("exportTables", self)

}


#' #' @export
#' exportTables.TrainTestResult <- function(self,
#'                                         export.dir = self$parameters$dir$tsv,
#'                                         #file.prefix,
#'                                          extension = ".tsv"){
#'
#'   message("\tExporting TrainTestResult object ", self$ID, " to tables")
#'   message("\t\tExport directory\t", export.dir)
#'
#'   dir.create(export.dir, showWarnings = FALSE, recursive = TRUE)
#'
#'  ## message("\t\tFile prefix\t", file.prefix)
#'   # UseMethod("exportTables", self)
#'
#'   ############## Exporting the count table ####################
#'   trainTestResult.file <- file.path(export.dir, paste( self$ID, "_Train_Test_Result", extension, sep = ""))
#'   message("\t\tExporting count table to TSV file\t",trainTestResult.file)
#'   write.table(self$stats, file = trainTestResult.file, row.names = FALSE, quote = FALSE, sep = "\t")
#'
#' NextMethod("exportTables", object)
#'
#'
#' }
