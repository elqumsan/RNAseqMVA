#' @title Create a StudyCase object, load RecountID dataset, and generate datasets.
#' @author Jacques van Helden and Mustafa AbuElQumsan
#' @description Create a StudyCase object, load RecountID dataset, and run preprocessing to generate the DataTableWitthClasses objects.

#' @param recountID A valide ID for an object of the ReCount2 database
#' @param parameters recountID-specific parameters specified in a YAML-formatted configuration file
#'
#' @return an object of class StudyCase, that is contain raw data before per-processing step, datasets for testing that are per-processed dataset and their respective parameters.
#'
#' @export

StudyCase  <- function (recountID, parameters) {


  message.with.time("Loading count table from recount", "; recountID = ", parameters$recountID)

  ## Run loadCounts method to get the counts per run + counts per sample + filtered counts
  result <- loadCounts(recountID = recountID,
                       parameters = parameters)

  ## Select training and testing sets on the filtered table with raw counts
  ## These wil then be passed to all the derived count tables (normalised, DGE, ...)

   result$filtered <- DataTableWithTrainTestSets(result$filtered)
   #result$filtered <-  UseMethod("DataTableWithTrainTestSets",result$filtered)

  ##### Normalize the counts without log2 transformation (first test) #####
  ##
  ## Note: this method takes a table with one column per sample and one
  ## row per gene, we thus have to transpose the raw count table.

  ###### Normalization method for the recount Table after merge and filtered it ########
  # dim(result$studyCases$norm$counts)
  message.with.time("Normalizing counts based on 75th percentile")
  result$norm <- NormalizeSamples(
    self = result$filtered,
    classColumn = result$parameters$classColumn,
    classColors = result$parameters$classColors,
    method = "quantile", quantile=0.75, log2 = FALSE)


  ##### Normalize counts with log2 transformation (second test) #####
  ##
  ## Note: this method takes a table with one column per sample and one
  ## row per gene, we thus have to transpose the raw count table.
  # if (parameters$compute) {
  message.with.time("Normalizing counts based on 75th percentile + log2 transformation")
  result$log2norm <- NormalizeSamples(
    self = result$filtered,
    classColumn = result$parameters$classColumn,
    method = "quantile", quantile=0.75,
    log2 = TRUE, epsilon=0.1)

  #### Derive an object having as features the principal components of log2norm ####

  ## Clone the log2norm object to copy all its parameters
  result$log2normPCs <- result$log2norm
  result$log2normPCs$dataType <- "log2normPCs"
  # names(result$log2normPCs)

  ## COmpte principal components
  result$log2normPCs$prcomp <-
    prcomp( t(na.omit(result$log2norm$dataTable)),
            center = TRUE,
            scale. = FALSE)
  ## Replace log2 normalised counts by principal components
  result$log2normPCs$dataTable <- t(result$log2normPCs$prcomp$x)
  # dim(result$log2norm$dataTable)
  # rownames(result$log2norm$dataTable)
  # dim(result$log2normPCs$dataTable)
  # rownames(result$log2normPCs$dataTable)
  # View(result$log2normPCs$dataTable)
  # biplot(result$log2normPCs$prcomp,cex=0.2) ## This is too heavy

  message.with.time("The computation of the variables importance by Random forest, and ordered it by the most importance")
  ## Clone the log2norm object to copy all its parameters
  result$log2normViRf <- result$log2norm
  result$log2normViRf$dataType <- "log2normViRf"

  rf.model  <- randomForest(
    x = t(result$log2norm$dataTable),
    y =  as.factor( result$log2norm$classLabels),
    xtest = t(result$log2norm$dataTable), importance = T, keep.forest = T)

  variable.importance <- importance(rf.model,type = 1,scale = F)
  #variable.importance[,1]

  ordered.varaible.importance <-order(variable.importance[,1],decreasing = T)

  ordered.dataTable.by.importace <-result$log2norm$dataTable[ordered.varaible.importance, ]

  sig.variables <- round(ncol(ordered.dataTable.by.importace) * 0.75)
  ordered.dataTable.by.importance  <- ordered.dataTable.by.importace[1:sig.variables, ]

  result$log2normViRf$viRf <- rf.model
  result$log2normViRf$ordereviRf <- ordered.varaible.importance
 # result$log2normViRf$sigviRf <- sig.variables
  result$log2normViRf$sigviRf <- ordered.dataTable.by.importance

  result$log2normViRf <- ordered.dataTable.by.importace

  ## Build a first version of the object based on passed parameters
  object <- structure(
    list(
      ID = recountID,
      parameters = parameters,
      rawData = list(
        countsPerRun = result$countsPerRun,
        countsPerSample = result$originalCounts),
      datasetsForTest = list(
        filtered = result$filtered,
        norm = result$norm,
        log2norm = result$log2norm,
        log2normPCs = result$log2normPCs,
        log2normViRf = result$log2normViRf
      )
    ),
    class="StudyCase")


  message("\t\tInstantiated an object of class StudyCase for recountID\t", recountID)
  return(object)
}

#' @title print a summary of an object belonging to class StudyCase
#' @author Mustafa AubElQumsan and Jacques van Helden
#' @description just print the summary of the object that is belonge to class StudyCase
#' @return print the summary of such object by utilizing generic function
#' @export
summary.StudyCase<- function(object){
  cat("\t\tObject belonge to StudyCase class\n")
  cat("\tRecountID        \t", object$ID ,"\n")
  cat("\t rawData = list( \n"  )
  cat("\t\tcountsPerRun = result$countsPerRun,\n")
  cat("\t\tcountsPerSample = result$originalCounts),\n")
  cat("\tdatasetsForTest = list( \n")
  cat("\t\tfiltered = result$filtered,\n")
  cat("\t\tnorm = result$norm,\n")
  cat("\t\tlog2norm = result$log2norm,\n")
  cat("\t\tlog2normPCs = result$log2normPCs)\n")
  cat("\t\tlog2normViRf = result$log2normViRf")
}

#' @export
summary.default<- function(object){
  cat("")
}

#' @title print a summary of an object belonging to class StudyCase
#' @author Mustafa AubElQumsan and Jacques van Helden
#' @description just print the summary of the object that is belonge to class StudyCase
#' @return print the summary of such object by utilizing generic function
#' @export
print.StudyCase <- function(object){
  summary(object)
}
