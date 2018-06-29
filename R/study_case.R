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
        log2normPCs = result$log2normPCs
      )
    ),
    class="StudyCase")


  message("\t\tInstantiated an object of class StudyCase for recountID\t", recountID)
  return(object)
}


summary.StudyCase<- function (object){

  cat("\t\tObject of", object$ID," belonges to StudyCases class "  , "\n")
  #cat("\tRecountID        \t", object$ID ,"\n")
  print( paste("rawDataTables are belonge to ", studyCases[[recountID]]$rawData$countsPerSample, "DataTablesWithClasses"))
}

<<<<<<< HEAD
=======


#' @title run variable importance by random Forest to on an object of class StudyCase.
#' @description run variable importance by random Forest to test importance on each feature of a data table.
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @param self object belong to StudyCase class.
#' @return an object of the same class as the input object
#' @export
RunViRf <- function(self) {
  message("\tRunning variable importance from Random Forest (RunViRf) for object of class ", paste(collapse=", ", class(self)))
  self <- UseMethod("RunViRf", self)
  return(self)
}


#' @title run variable importance by random Forest on an object of class StudyCase
#' @description run variable importance by random Forest on an object of class StudyCase to test importance on each feature of a data table, and order variables by decreasing the importance of features in dataTable.
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @param self object belong to StudyCase class.
#' @return a clone of the input StudyCase object with an added
#' DataTableWithTrainTestSets containing the log2norm data table
#' where features have been re-ordered by decreasing the importance of the features.
#' @export
RunViRf.StudyCase <- function(self) {
  message.with.time("Defining gene order according to variable importance by random Forest")

  self$datasetsForTest$log2norm_ViRf_sorted <- self$datasetsForTest$log2norm
  self$datasetsForTest$log2norm_ViRf_sorted$dataType <- "log2normViRf"

  rf.model  <- randomForest(
    x = t(self$datasetsForTest$log2norm$dataTable),
    y =  as.factor( self$datasetsForTest$log2norm$classLabels),
    xtest = t(self$datasetsForTest$log2norm$dataTable), importance = T, keep.forest = T)
  variable.importance <- importance(rf.model, type = 1, scale = F)
  ordered.varaible.importance <-order(variable.importance[,1],decreasing = T)

  ordered.dataTable.by.importace <-self$datasetsForTest$log2norm$dataTable[ordered.varaible.importance, ]
  sig.variables <- round(nrow(ordered.dataTable.by.importace) * 0.75)
  ordered.dataTable.by.importance  <- ordered.dataTable.by.importace[1:sig.variables, ]

  self$datasetsForTest$log2norm_ViRf_sorted$viRf <- rf.model
  self$datasetsForTest$log2norm_ViRf_sorted$ordereviRf <- ordered.varaible.importance
  self$datasetsForTest$log2norm_ViRf_sorted$sigviRf <- ordered.dataTable.by.importance
  self$datasetsForTest$log2norm_ViRf_sorted$orderedDataTable <- ordered.dataTable.by.importace
  self$datasetsForTest$log2norm_ViRf_sorted$dataTable <- ordered.dataTable.by.importace
  return(self)

}


#' @title Draw an histogram with the raw  counts per gene
#' @author Jacques van Helden & Mustafa AbuElQumsan
histCountsPerGeneClass.StudyCase <- function(self) {
  counts <- as.vector(unlist(self$datasetsForTest$filtered$dataTable))
  hist(log2(counts + epsilon), breaks=100,
       main = paste(recountID,
                    " – Histogram of log2(counts)"),
       xlab = "log2(counts)",
       ylab = "Number of genes",
       col="#CCBBFF")
  legend("topright",
         legend = paste("Max counts per gene =",
                        prettyNum(max(counts), big.mark = ",")))
}

#' @title Draw an XY plot to compare the mean counts per gene for two user-specified class
XYplot.StudyCase <- function(self,
                             class1,
                             class2) {

  normcounts <- self$datasetsForTest$norm$dataTable
  gene.mean.per.class <- by(
    t(normcounts),
    INDICES = self$datasetsForTest$norm$classLabels, FUN = colMeans)

  epsilon <- 0.1
  x1 <- gene.mean.per.class[[class1]] + epsilon
  x2 <- gene.mean.per.class[[class2]] + epsilon


  ## XY plot
  plot (x = log2(x1),
        y = log2(x2),
        main = paste(recountID, "\nlog2(scaled counts) per gene"),
        xlab = "Bone marrow ",
        ylab = "Heparinised blood",
        las = 1,
        col = densCols(x = log2(x1), y = log2(x2)),
        panel.first = grid())
  abline(a = 0, b = 1, col = "black", lwd = 1)

  ## MA plot
  A <- (log2(x1) + log2(x2))/2
  M <- log2(x1) - log2(x2)
  plot (x = A, y = M,
        main = paste(recountID, "\nMA plot"),
        xlab = "A", ylab = "M",
        col = densCols(x = A, y = M),
        panel.first = grid())
  abline(h = 0, col = "black", lwd = 1)

}

>>>>>>> f2616f4f21f0586c16190bc4b2ba7826ad19420d
