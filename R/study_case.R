#' @title Create a StudyCase object, load RecountID dataset, and generate datasets.
#' @author Jacques van Helden and Mustafa AbuElQumsan
#' @description Create a StudyCase object, load RecountID dataset, and run preprocessing to generate the DataTableWitthClasses objects.
#' @param recountID a valid experiment ID for a record of the ReCount2 database
#' @param parameters recountID-specific parameters specified in a YAML-formatted configuration file
#'
#' @return an object of class StudyCase, which contains raw data (before preprocessing steps),
#' as well as several datasets normalised with different methods.
#'
#' @export

StudyCase  <- function(recountID, parameters) {

  if (is.null(parameters$standardization$epsilon)) {
    epsilon <- 0.1
    message("epsilon (pseudo-count) not defined in config file -> using default value (", epsilon, ")")
  } else {
    epsilon <- parameters$standardization$epsilon
  }

  message.with.time("Loading count table from recount", "; recountID = ", parameters$recountID)

  ## Run loadCounts method to get the counts per run + counts per sample + filtered counts
  result <- loadCounts(recountID = recountID, parameters = parameters)
  # class(result$filtered)

  ## Cast the filtered counts from class DataTableWithClasses to class DataTableWithTrainTestSets.
  ## This will also append attributes with selected training and testing sets on the filtered table with raw counts
  ## These wil then be passed to all the derived count tables (normalised, DGE, ...).
  result$filtered <- DataTableWithTrainTestSets(result$filtered)

  ##### Normalize the counts without log2 transformation (first test) #####
  ##
  ## Note: this method takes a table with one column per sample and one
  ## row per gene, we thus have to transpose the raw count table.

  ## Get standardization methods and options from parameters
  standardization.methods <- parameters$standardization$method
  quantile <- parameters$standardization$quantile
  dataNames <- vector()
  for (method in standardization.methods) {

    ###### Normalization method for the recount Table after merge and filtered it ########
    # dim(result$studyCases$norm$counts)
    if (method == "quantile") {
      method.name <- paste(sep = "", "q", quantile)
    } else {
      method.name <- method
    }
    result[[method.name]] <- NormalizeSamples(
      self = result$filtered,
      method = method,
      quantile = quantile,
      log2 = FALSE)
    dataNames <- append(dataNames, method.name)

    ##### Normalize counts by log2 transformation #####
    message("\t\tComputing log2-transformed normalised counts")
    log2.name <- paste(sep = "", method.name, "_log2")
    result[[log2.name]] <- result[[method.name]]
    result[[log2.name]][["dataType"]] <- log2.name
    result[[log2.name]][["dataTable"]] <- log2(result[[method.name]][["dataTable"]] + epsilon)
    dataNames <- append(dataNames, log2.name)

    #### Derive an object having as features the principal components of log2norm ####
    message("\t\tComputing principal components")
    pc.name <- paste(sep = "", log2.name, "_PC")
    result[[pc.name]] <- result[[log2.name]] ## Clone the log2norm object to copy all its parameters
    result[[pc.name]][["dataType"]] <- pc.name
    dataNames <- append(dataNames, pc.name)

    counts.for.pc <- result[[pc.name]]$dataTable # Suppress rows containing NA values

    ## Filter out samples with Infinite values
    ## (results from null size factor -> infinite scaling factor)
    samples.with.inf <- apply(is.infinite(counts.for.pc), 2, sum) > 0
    if (sum(samples.with.inf) > 0) {
      message("\tDiscarding ", sum(samples.with.inf), " samples with infinite values (due to null size factor)")
      counts.for.pc <- counts.for.pc[, !samples.with.inf]
      message("\t\tRemaining: ",
              nrow(counts.for.pc), " features (rows) x ",
              ncol(counts.for.pc), " samples (columns).")


    }

    ## Remove features with NA values
    features.with.NA <- apply(is.na(counts.for.pc), 1, sum) > 0
    if (sum(features.with.NA) > 0) {
      message("Discarding ", sum(features.with.NA), " samples with infinite values (due to null size factor)")
      counts.for.pc <- na.omit(counts.for.pc) # Suppress rows containing NA values
      message("\t\tRemaining: ",
              nrow(counts.for.pc), " features (rows) x ",
              ncol(counts.for.pc), " samples (columns).")
    }

    ## Compute principal components
    result[[pc.name]]$prcomp <-
      prcomp(t(counts.for.pc),
             center = TRUE,
             scale. = FALSE)
    ## Replace log2 normalised counts by principal components
    result[[pc.name]]$dataTable <- t(result[[pc.name]]$prcomp$x)
    # names(result[[pc.name]])
    # hist(result[[pc.name]][["dataTable"]], breaks=200)

  }

  ## Build a first version of the object based on passed parameters
  self <- structure(
    list(
      ID = recountID,
      parameters = parameters,
      rawData = list(
        countsPerRun = result$countsPerRun,
        countsPerSample = result$originalCounts),
      datasetsForTest = list(
        filtered = result$filtered)
    ),
    class = "StudyCase")

  ## Append the normalised data types to the StudyCase object
  for (name in dataNames) {
    self$datasetsForTest[[name]] <- result[[name]]
  }
  self$normalisedDataNames <- dataNames

  #### Define variable ordering based on DESeq2 significance ####
  if ("DESeq2" %in% project.parameters$global$ordering.methods) {
    self <- RunDESeq2(self)
  }

  #### Define variable (gene) ordering according to edgeR significance ####
  if ("edgeR" %in% project.parameters$global$ordering.methods) {
    message.with.time("Defining gene order according to edgeR differential expression")
    self <- RunedgeR(self)
  }


  #### Define variable (gene) ordering according to the variable importanve (VI) computed by RandomForest ####
  if ("RF" %in% project.parameters$global$ordering.methods) {
    message.with.time("Computing variables importance by Random Forest (RF), and ordering features by decreasing importance. ")
    # ## Clone the log2norm object to copy all its parameters
    self <- RunViRf(self)
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
#' @description run edgeR on an object of class StudyCase to test differential expression on each feature of a data table, and order variables by increasing adjusted p-value.
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @param self object belong to StudyCase class
#' @return an object of the same class as the input object
#' @export
RunedgeR <- function(self) {
  message("\tRunning edgeR for object of class ", paste(collapse=", ", class(self)))
  self <- UseMethod("RunedgeR", self)
  return(self)
}


#' @title run edgeR on an object of class StudyCase
#' @description run edgeR on an object of class StudyCase to test differential expression on each feature of a data table, and order variables by increasing adjusted p-value.
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @param self object belong to StudyCase class
#' @return a clone of the input StudyCase object with an added
#' DataTableWithTrainTestSets containing the log2norm data table
#' where features have been re-ordered by increasing adjusted p-value.
#' @export
RunedgeR.StudyCase <- function(self) {

  message.with.time("Defining gene order according to edgeR differential expression")

  self$datasetsForTest$log2norm_edgeR_sorted <- self$datasetsForTest$log2norm

  ## include the edgeR result table in the resulting data object
  self$datasetsForTest$log2norm_edgeR_sorted$edgeR  <-
    DEGordering(countDataset = self$datasetsForTest$filtered, method = "edgeR")

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
#' @description run DESeq2 on an object of class StudyCase to test differential expression on each feature of a data table, and order variables by increasing adjusted p-value.
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @param self object belong to StudyCase class
#' @return an object of the same class as the input object
#' @export
RunDESeq2 <- function(self) {
  message("\tRunning DESeq2 for object of class ", paste(collapse=", ", class(self)))
  self <- UseMethod("RunDESeq2", self)
  return(self)
}


#' @title run DESeq2 on an object of class StudyCase
#' @description run DESeq2 on an object of class StudyCase to test differential expression on each feature of a data table, and order variables by increasing adjusted p-value.
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @param self object belong to StudyCase class.
#' @return a clone of the input StudyCase object with an added
#' DataTableWithTrainTestSets containing the log2norm data table
#' where features have been re-ordered by increasing adjusted p-value.
#' @export
RunDESeq2.StudyCase <- function(self) {

  message.with.time("Defining gene order according to DESeq2 differential expression")

  self$datasetsForTest$log2norm_DESeq2_sorted <- self$datasetsForTest$log2norm

  ## include the DESeq2 result table in the resulting data object
  self$datasetsForTest$log2norm_DESeq2_sorted$DESeq2  <-
    DEGordering(self$datasetsForTest$filtered, method = "DESeq2")

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



#' @title Compute variable importance with Random Forest on an object of class StudyCase, and
#' @description run variable importance by random Forest to test importance on each feature of a data table.
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @param self an object belonging to StudyCase class.
#' @return a clone of the input StudyCase object with an added
#' DataTableWithTrainTestSets containing the log2norm data table
#' where features have been re-ordered by decreasing the importance of the features.
#' @export
RunViRf <- function(self) {
  message("\tRunning variable importance from Random Forest (RunViRf) for object of class ", paste(collapse=", ", class(self)))
  self <- UseMethod("RunViRf", self)
  return(self)
}


#' @title Sort variables by decreasing variable importance as computed by Random Forest
#' @description run Random Forest on an object of class StudyCase, and, and order a dataTable by decreasing order of RF variable importance.
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @param self an object belonging to StudyCase class.
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


#' @title Draw an histogram with the raw counts per gene
#' @description to carify the distribution of the raw counts per gene.
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @param self object belong to StudyCase class.
#' @return an object of the same class as the input object
#' @export
histCountsPerGeneClass <- function(self) {
  message("\tRunning Drawing an histogram with the raw count per gene for object of class ", paste(collapse=", ", class(self)))
  # self <- UseMethod("histCountsPerGeneClass", self)
  UseMethod("histCountsPerGeneClass", self)
  # return(self)
}

#' @title Draw an histogram with the raw  counts per gene
#' @author Jacques van Helden & Mustafa AbuElQumsan
#' @param self object belong to StudyCase class.
#' @return an object of the same class as the input object
#' @export
histCountsPerGeneClass.StudyCase <- function(self) {
  counts <- as.vector(unlist(self$datasetsForTest$filtered$dataTable))
  hist(log2(counts + epsilon), breaks=100,
       main = paste(recountID,
                    " - Histogram of log2(counts)"),
       xlab = "log2(counts)",
       ylab = "Number of genes",
       col="#CCBBFF")
  legend("topright",
         legend = paste("Max counts per gene =",
                        prettyNum(max(counts), big.mark = ",")))
}

#' @title Draw an XY plot to compare the mean counts per gene for two user-specified class
#' @description Drawing an XY plot to compare the mean counts per gene for two user-specified class of the raw counts per gene.
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @param self object belong to StudyCase class.
#' @return an object of the same class as the input object
#' @export
XYplot <- function(self) {
  message("\tRunning Drawing an XY plot to compare the mean counts per gene with the raw count per gene for object of class ", paste(collapse=", ", class(self)))
  # self <- UseMethod("histCountsPerGeneClass", self)
  UseMethod("XYplot", self)
  # return(self)
}

#' @title Draw an XY plot to compare the mean counts per gene for two user-specified class
#' @description Drawing an XY plot to compare the mean counts per gene for two user-specified class of the raw counts per gene.
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @param self object belong to StudyCase class.
#' @param class1 the name of the first class to plot
#' @param class2 the name of the second class to plot
#' @return an object of the same class as the input object
#' @export
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

