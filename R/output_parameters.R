#' @title Compute output parameters for an RNAseqMVA analysis
#' @description Compute output parameters (output directory, file prefix, figure label)  for a given analysis given its input parameterrs (recountIR, dataset, classifier, permutation or not).
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @param dataset an object of class DataTableWithTrainTestSets
#' @param classifier is the type of classifier that is used with repeated process.
#' @param permute=FALSE is show if the class lable are permuted this for sake of the knowing the strength and weaknesses of the classifier
#' @param createDir=TRUE if TRUE, the result directory is created automatically if it does not exist
#'
#' @return a list with the following fields
#' \itemize{
#'   \item resultDir: path to the result directory, in which the result files should be stored.
#'   \item filePrefix: string that should prepend the name of each result file.
#'   \item figureLabel:  short string with a label for the graphics
#' }
#' @export
outputParameters <- function(dataset,
                             classifier,
                             permute = FALSE,
                             createDir = TRUE,
                             knn.k = NULL) {

  ## Check the class of hte dataset
  if (!is(object = dataset, class2 = "DataTableWithTrainTestSets")) {
    stop("filePrefix requires an object of class DataTableWithTrainTestSets")
  }

  #  Verbosity
  message("\tDefining file Prefix for dataset ", dataset$ID,
          "; data type: ", dataset$dataType,
          "; classifier: ", classifier)

  ## Define the piece of prefix and output dir that indicates the classifier

  ## For some classifiers, the file prefix includes the main parameters
  if (classifier == "svm") {
    ## For SVM, the prefix includes the kernel
    if (is.null(dataset$parameters$svm$kernel)) {
      dataset$parameters$svm$kernel = "linear"
    }
    classifierPrefix <- paste(sep = "", "svm_", dataset$parameters$svm$kernel)

  } else if (classifier == "knn") {
    ## For KNN, the prefix includes the value of the k parameter
    if (is.null(knn.k)) {
      knn.k <- dataset$parameters$knn$k
    }
    classifierPrefix <- paste(sep = "", "knn_k", knn.k)

  } else {
    ## For all other types of calssifer, the variable
    ## classifierPrefix is simply the name of the classifier
    classifierPrefix <- classifier
  }


  ## Define the path to the result directory
  resultDir <- file.path(dataset$parameters$dir$results,
#                         paste(sep = "_", dataset$ID, dataset$parameters$feature),
                         classifier)
  resultDir <- gsub(pattern = " ", replacement = "_", x = resultDir) ## Avoid spaces in file names

  if (createDir) {
    dir.create(resultDir, showWarnings = FALSE, recursive = TRUE)
  }

  ## Define the path to the tables TSV result directory
  resultDirTablesTSV <- file.path(dataset$parameters$dir$results,
                         # dataset$ID,
                         classifier, "tables")
  resultDirTablesTSV <- gsub(pattern = " ", replacement = "_", x = resultDirTablesTSV) ## Avoid spaces in file names
  dir.create(resultDirTablesTSV, showWarnings = FALSE, recursive = TRUE)

  ## Define file prefix
  filePrefix <- paste(sep = "_", dataset$ID,
                      dataset$parameters$feature,
                      classifierPrefix,
                      dataset$dataType)
  if (permute) {
    filePrefix <- paste(sep = "_", filePrefix, dataset$parameters$perm.prefix)
  }
  filePrefix <- gsub(pattern = " ", replacement = "_", x = filePrefix) ## Avoid spaces in file names

  ## Define file label for figures
  figLabel <- paste(sep = " ", classifierPrefix,  dataset$dataType)
  if (permute) {
    figLabel <- paste(sep = " ", figLabel, dataset$parameters$perm.prefix)
  }
  figLabel <- gsub(pattern = "_", replacement = " ", x = figLabel) ## For figure labels spaces are more readable than underscores

  # message("\t\tfilePrefix\t", filePrefix)
  # message("\t\tresultDir\t", resultDir)
  # message("\t\tresultDirTablesTSV\t", resultDirTablesTSV)
  # message("\t\tfigLabel\t", figLabel)


  return(list("resultDir" = resultDir,
              "resultDirTablesTSV" = resultDirTablesTSV,
              "filePrefix" = filePrefix,
              "figLabel" = figLabel))
}
