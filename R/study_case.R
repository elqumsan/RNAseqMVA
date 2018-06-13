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

  ## Compute principal components
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


  ##### instantiate object from ged-dataSet from Differential analysis with DESeq2 and edgeR to define gene (variable) order ####
  # message.with.time("instantiate object of Differential analysis with DESeq2 and edgeR to define gene (variable) order")
  # result$DEGdataSets <- result$filtered
  # #result$DEGdataSets$edgeR <- list()
  # result$DEGdataSets$DESeq2  <- DEGordering(result$originalCounts$dataTable, result$originalCounts$classLabels,
  #                                           method = project.parameters$global$ordering.methods[1] , randomized = TRUE )
  # result$DEGdataSets$DESeq2$dataType <- "DESeq2orderedVariables"
  #
  # result$DEGdataSets$edgeR  <- DEGordering(result$originalCounts$dataTable, result$originalCounts$classLabels,
  #                                          method = project.parameters$global$ordering.methods[2] , randomized = TRUE )
  # result$DEGdataSets$edgeR$dataType <- "edgeRorderedVariables"

   ## Build a first version of the object based on passed parameters
  self <- structure(
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
        log2normPCs = result$log2normPCs
#        log2normViRf = result$log2normViRf# ,
#        log2norm_DESeq2_sorted = result$log2norm_DESeq2_sorted,
#        log2norm_edgeR_sorted = result$log2norm_edgeR_sorted
      )
    ),
    class="StudyCase")

  #### DESeq2-sorted variables ####
  if ("DESeq2" %in% project.parameters$global$ordering.methods) {
    RunDESeq2(self)
    # message.with.time("Defining gene order according to DESeq2 differential expression")
    #
    # self$log2norm_DESeq2_sorted <- self$datasetsForTest$log2norm
    # self$log2norm_DESeq2_sorted$DESeq2  <-
    #   DEGordering(self$datasetsForTest$filtered,
    #               method = "DESeq2", randomized = TRUE )
    # self$DEGdataSets$DESeq2$dataType <- "log2norm_DESeq2_ordered"
    #
    # ## Note: we use the log2norm as variables,
    # ## but sort them according to DESeq2,
    # ## which was based on the raw counts
    # self$log2norm_DESeq2_sorted$dataTable <-
    #   self$log2norm_DESeq2_sorted$dataTable[self$log2norm_DESeq2_sorted$DESeq2$geneOrder, ]
  }

  #### edgeR-sorted variables ####
  if ("edgeR" %in% project.parameters$global$ordering.methods) {
    message.with.time("Defining gene order according to edgeR differential expression")

    self$log2norm_edgeR_sorted <- self$datasetsForTest$log2norm
    self$log2norm_edgeR_sorted$edgeR  <-
      DEGordering(self$datasetsForTest$filtered,
                  method = "edgeR", randomized = TRUE )
    self$DEGdataSets$edgeR$dataType <- "log2norm_edgeR_ordered"

    ## Note: we use the log2norm as variables,
    ## but sort them according to edgeR,
    ## which was based on the raw counts
    self$log2norm_edgeR_sorted$dataTable <-
      self$log2norm_edgeR_sorted$dataTable[self$log2norm_edgeR_sorted$edgeR$geneOrder, ]
  }


  if ("RF" %in% project.parameters$global$ordering.methods) {

    message.with.time("Computing variables importance by Random Forest (RF), and ordering features by decreasing importance. ")
    # ## Clone the log2norm object to copy all its parameters


    self$log2norm_ViRf_sorted <- self$datasetsForTest$log2norm
    self$log2norm_ViRf_sorted$dataType <- "log2normViRf"
    rf.model  <- randomForest(
      x = t(result$log2norm$dataTable),
      y =  as.factor( result$log2norm$classLabels),
      xtest = t(result$log2norm$dataTable), importance = T, keep.forest = T)
    variable.importance <- importance(rf.model, type = 1, scale = F)
    ordered.varaible.importance <-order(variable.importance[,1],decreasing = T)
    ordered.dataTable.by.importace <-result$log2norm$dataTable[ordered.varaible.importance, ]
    sig.variables <- round(nrow(ordered.dataTable.by.importace) * 0.75)
    ordered.dataTable.by.importance  <- ordered.dataTable.by.importace[1:sig.variables, ]
    self$log2norm_ViRf_sorted$viRf <- rf.model
    self$log2norm_ViRf_sorted$ordereviRf <- ordered.varaible.importance
    self$log2norm_ViRf_sorted$sigviRf <- ordered.dataTable.by.importance
    self$log2norm_ViRf_sorted$orderedDataTable <- ordered.dataTable.by.importace
    self$log2norm_ViRf_sorted$dataTable <- ordered.dataTable.by.importace
  }

  message("\t\tInstantiated an object of class StudyCase for recountID\t", recountID)
  return(self)
}

#' @title print a summary of an object belonging to class StudyCase
#' @author Mustafa AubElQumsan and Jacques van Helden
#' @description just print the summary of the object that is belonge to class StudyCase
#' @return print the summary of such object by utilizing generic function
#' @export
summary.StudyCase<- function(self){
  cat("StudyCase object\n")
  cat("\tRecountID\t", self$ID ,"\n")
  cat("\trawData\n")
  for (dataset in names(self$rawData)) {
    cat ("\t\t", dataset, "\t", paste(collapse=", ", class(self$rawData[[dataset]])), "\n")
  }
  cat("\tdatasetsForTest\n")
  for (dataset in names(self$datasetsForTest)) {
    cat ("\t\t", dataset, "\t", paste(collapse=", ", class(self$datasetsForTest[[dataset]])), "\n")
  }
}

#' @export
summary.default<- function(self){
  cat("")
}

#' @title print a summary of an object belonging to class StudyCase
#' @author Mustafa AubElQumsan and Jacques van Helden
#' @description just print the summary of the object that is belonge to class StudyCase
#' @return print the summary of such object by utilizing generic function
#' @export
print.StudyCase <- function(self){
  summary(self)
}



#' @title run edgeR to test differential expression on each feature of a data table.
#' @description
#' @param self object
#' @return an object of the same class as the input object
#' @export
RunedgeR <- function(self) {
  message("\tRunning edgeR for object of class ", paste(collapse=", ", class(self)))
  UseMethod("RunedgeR", self)
  return(self)
}


#' @title run edgeR on an object of class StudyCase
#' @description run edgeR on an object of class StudyCase to test differential expression on each feature of a data table, and order variables by increasing adjusted p-value.
#' @param self object
#' @return a clone of the input StudyCase object with an added
#' DataTableWithTrainTestSets containing the log2norm data table
#' where features have been re-ordered by increasing adjusted p-value.
#' @export
RunedgeR.StudyCase <- function(self) {

  message.with.time("Defining gene order according to edgeR differential expression")

  self$datasetsForTest$log2norm_edgeR_sorted <- self$datasetsForTest$log2norm

  ## include the edgeR result table in the resulting data object
  self$datasetsForTest$log2norm_edgeR_sorted$edgeR  <-
    DEGordering(self$datasetsForTest$filtered,
                method = "edgeR", randomized = TRUE )

  ## Specify te data type
  self$datasetsForTest$log2norm_edgeR_sorted$dataType <- "log2norm_edgeR_ordered"

  ## Note: we use the log2norm as variables,
  ## but sort them according to edgeR,
  ## which was based on the raw counts
  self$datasetsForTest$log2norm_edgeR_sorted$dataTable <-
    self$datasetsForTest$log2norm_edgeR_sorted$dataTable[
      self$datasetsForTest$log2norm_edgeR_sorted$edgeR$geneOrder, ]

  return(self)
}


#' @title run DESeq2 to test differential expression on each feature of a data table.
#' @description
#' @param self object
#' @return an object of the same class as the input object
#' @export
RunDESeq2 <- function(self) {
  message("\tRunning DESeq2 for object of class ", paste(collapse=", ", class(self)))
  UseMethod("RunDESeq2", self)
  return(self)
}


#' @title run DESeq2 on an object of class StudyCase
#' @description run DESeq2 on an object of class StudyCase to test differential expression on each feature of a data table, and order variables by increasing adjusted p-value.
#' @param self object
#' @return a clone of the input StudyCase object with an added
#' DataTableWithTrainTestSets containing the log2norm data table
#' where features have been re-ordered by increasing adjusted p-value.
#' @export
RunDESeq2.StudyCase <- function(self) {

  message.with.time("Defining gene order according to DESeq2 differential expression")

  self$datasetsForTest$log2norm_DESeq2_sorted <- self$datasetsForTest$log2norm

  ## include the DESeq2 result table in the resulting data object
  self$datasetsForTest$log2norm_DESeq2_sorted$DESeq2  <-
    DEGordering(self$datasetsForTest$filtered,
                method = "DESeq2", randomized = TRUE )

  ## Specify te data type
  self$datasetsForTest$log2norm_DESeq2_sorted$dataType <- "log2norm_DESeq2_ordered"

  ## Note: we use the log2norm as variables,
  ## but sort them according to DESeq2,
  ## which was based on the raw counts
  self$datasetsForTest$log2norm_DESeq2_sorted$dataTable <-
    self$datasetsForTest$log2norm_DESeq2_sorted$dataTable[
      self$datasetsForTest$log2norm_DESeq2_sorted$DESeq2$geneOrder, ]

  return(self)
}



