#' @title return output parameters (output directory, file prefix, figure label)  for a given analysis given its parameterrs (dataset, classifier, permutation or not)
#' @description for saking of saveing all required output parameters (output dirctory, file prefix and figure label) to facilitate and reduce the  process of plotting
#' and saving the results in identified workspace.
#' @author Mustafa ABUELQUMSAN and Jacques van Helden
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
                             createDir = TRUE) {

  # ## Check the class of hte dataset
  # if (!is(object = dataset, class2 = "DataTableWithTrainTestSets")) {
  #   stop("filePrefix requires an object of class DataTableWithTrainTestSets")
  # }

  ##  Verbosity
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
    classifier_prefix = paste(sep = "", "svm_", dataset$parameters$svm$kernel)
  } else if (classifier == "knn") {
    ## For KNN, the prefix includes the value of the k parameter
    if (is.null(dataset$parameters$knn$k)) {
      dataset$parameters$knn$k
    }
    classifier_prefix = paste(sep = "", "knn_k", dataset$parameters$knn$k)

  } else {
    ## For all other types of calssifer, the variable
    ## classifier_prefix is simply the name of the classifier
    classifier_prefix = classifier
  }


  ## Define the path to the result directory
  resultDir <- file.path(dataset$parameters$dir$results,
                         dataset$ID,
                         classifier)

  dir.create(resultDir, showWarnings = FALSE, recursive = TRUE)

  ## Define the path to the tables TSV result directory
  resultDirTablesTSV <- file.path(dataset$parameters$dir$results,
                         dataset$ID,
                         classifier, "tables")
  dir.create(resultDirTablesTSV, showWarnings = FALSE, recursive = TRUE)

  ## Define file prefix
  filePrefix <- paste(sep = "_", dataset$ID, classifier_prefix,  dataset$dataType)
  if (permute) {
    filePrefix <- paste(sep = "_", filePrefix, dataset$parameters$perm.prefix)
  }
  filePrefix <- sub(pattern = " ", replacement = "_", x = filePrefix) ## Avoid spaces in file names

  ## Define file label for figures
  figLabel <- paste(sep = " ", classifier_prefix,  dataset$dataType)
  if (permute) {
    figLabel <- paste(sep = " ", figLabel, dataset$parameters$perm.prefix)
  }
  figLabel <- sub(pattern = "_", replacement = " ", x = figLabel) ## For figure labels spaces are more readable than underscores

  message("\t\tfilePrefix\t", filePrefix)


  return(list("resultDir" = resultDir,
              "resultsTablesTSV" =resultDirTablesTSV,
              "filePrefix" = filePrefix,
              "figLabel" = figLabel))
}
