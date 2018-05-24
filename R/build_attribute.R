#' @title build attributes for an object depending on its class
#' @author Mustafa AbuELQumsan and Jacques van Helden
#' @description building all the attributes for the object based on its own class. in order to give it some certian features based on its own class e.g. we have DataTableWithClases, dataTableWithTrainTestSets, ect.
#' @param self which is dataTable (a data.frame) with one row per feature (genes) and one column per individuals (sample)
#'
#' @return an object that has all required attributes as in the following attributes.
#' \itemize {
#'  \item nbSamples:  that is the number of samples in relatedness object.
#'  \item nbGenes:    that is the number of genes in  relatedness object.
#'  \item sampleNames: that is the names of all samples involved in related object.
#'  \item geneNames:   that is names of all genes in the relatedness object.
#'  \item classLabels: that is the class labels for all the samples in the relatedness object.
#'  \item classNames:  that is the names of all classes included in the relatedness object.
#'  \item nbClasses:   that is the overall number of the classes involved in the relatedness object.
#'  \item classProperities: that is the data.frame composed of two columns one for the class name, number of samples in each class, train size per class
#'  \item samplesPerClass: such paraperter are involved in classProperities.
#'  \item classFrequencies: such parameter means the ratio of certain samples for somehow class so as to the whole samples into the dataTable.
#'  \item randExpectedHitRate: it is corss product for the classFrequensies .
#'  \item randExpectedMisclassificationRate: it is 1 - randExpectedHitRate.
#'  \item sampleColors: such parameter for visualisation of the sample by putting some useful colors.

#' }
#'
#' @export
buildAttributes <- function(self) {
  message("\tBuilding attributes for object of class ", paste(collapse=", ", class(self)))
  UseMethod("buildAttributes", self)
#  message("trainIndices length : ", length(self$trainTestProperties$trainIndices))
#  message("\t\treturning from buildAttributes()")
  return(self) ## NOT SURE THIS LINE IS REQUIRED. TO BE CHECKED.
}


#' @title default method to build attributes for an object.
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @description Just send message with object classes. The class-specific builders should have been be called before.
#' @param self that is an object ahs been instantiated.
#' @return print message that inform the instantiate the object has been done.
#'
#' @export
buildAttributes.default <- function(self) {
  message("\tFinished building attributes for object of class ", paste(collapse=",", class(self)))
#  print (self$trainTestProperties)
#  message("trainIndices length : ", length(self$trainTestProperties$trainIndices))
#  message("\t\treturning from buildAttributes.default()")
  return(self)
}


