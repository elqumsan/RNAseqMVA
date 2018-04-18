#' @title file prefix procedures with a given classifier and a given object
#' @author Mustafa ABUELQUMSAN and Jacques van Helden
#' @description due to such file will be called from other scripts in order to computes some classifier-specific parameters
#' @param self an object of class DataTableWithTrainTestSets
#' @param classifier is the type of classifier that is used with repeated process.
#' @param permute is show if the class lable are permuted this for sake of the knowing the strength and weaknesses of the classifier
#' @param file.prefix in order to let us to save file for any other script
#'
#' @return
#' @import foreach
#' @import doParallel
#' @export


filePrefix <- function( self,
                        classifier,
                        permute = FALSE,
                        file.prefix = NULL){

  message("\tDefine the file Prefix for each classifier for saking study the impact of ",classifier," parameters:")
  ## TO DO:THIS SHOULD BE MOVED TO A SEPARATE FUNCTION
  ## Define file prefix is not specified in paramters
  if (is.null(file.prefix)) {
    if (classifier == "svm") {
      if (is.null(parameters$svm$kernel)) {
        parameters$svm$kernel = "linear"
      }
      classifier_prefix = paste(sep="", "svm_", self$parameters$svm$kernel)
    } else if (classifier == "knn") {
      if (is.null(parameters$knn$k)) {
        parameters$knn$k
      }
      classifier_prefix = paste(sep="", "knn_", self$parameters$knn$k)

    } else {
      classifier_prefix = classifier
    }
    file.prefix <- paste(sep="_", classifier_prefix, self$ID,  self$dataType, self$variablesType)
    file.prefix <- sub(pattern = " ", replacement = "_", x = file.prefix) ## Avoid spaces in file names
    if (permute) {
      file.prefix <- paste( file.prefix, "permLabels")
    }
  }

return(file.prefix)
}
