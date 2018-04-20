#' @title constructor of the DataTableWithClasses class
#' @author Mustafa AbuElQumsan and Jacques van helden
#' @description This class contains count tables
#' (typically used to represent expression profiles in RNA-seq analysis)
#' documented with a Pheno table (description of each sample) and completed
#' with some additional fields for convenience: sample names, gene names, ...
#' @param dataTable a data.frame with one row per feature (e.g. gene) and
#' one column per sample.
#' @param phenoTable a data frame describing each sample: one row per sample and one column per attribute of a sample.
#' @param dataType="raw counts" data type, free text (e.g. raw counts, log2-transformed counts, log2 normalised counts  ...).
#' @param parameters global and specific parameters for the analysis of this recountID
#'
#' @export

DataTableWithClasses <- function(dataTable,
                                  phenoTable,
                                  # classColumn = parameters$classColumn,
                                  # classColors = parameters$classColors,
                                  # ID = parameters$recountID,
                                  # variablesType,
                                  dataType,
                                  parameters) {

  ## Check required parameters
  for (p in c("recountID", "classColumn")) {
    if (is.null(parameters[[p]])) {
      stop("Missing required parameter: '", p,
           "'.\n\tPlease check configuration file. ")
    } else {
      assign(p, parameters[[p]])
    }
  }

  ## Built a list from the input parameters
  message.with.time("\tCreating object of class DataTableWithClasses\t", parameters$recountID)

  ## Build a first version of the object based on passed parameters
  object <- structure(
    list(
      ID = recountID,
      dataTable = dataTable,
      phenoTable = phenoTable,
      # variablesType = variablesType,
      dataType = dataType,
      parameters = parameters,

      ## NOTE: classColumn and classColors are now attached to the object via parameters
      ## -> should be suppressed from below but then we need to check that everything works fine.
      classColumn = classColumn# ,
#      classColors = classColors
    ),
    class="DataTableWithClasses")
  # names(object)
  # class(object)
  # attributes(object)
  # UseMethod("exportTables", self)

  ##  Get user-specificed class colors if defined. If not, will be defined with buildAttributes
  if (!is.null(parameters$classColors)) {
    object$classColors <- parameters$classColors
  }



  ## Define all the derived attributes from the count table and pheno table
   object <- buildAttributes(object)
  # object <- UseMethod("buildAttributes", object)
  ## Attach the parameters to the object
#  object$parameters <- parameters

  message("\t\tInstantiated an object of class dataTablewithClasses for recountID\t", recountID)
  return(object)
}

#' @title print a summary of an object belonging to class DataTableWithClasses
#' @author Mustafa AubElQumsan and Jacques van Helden
#' @export
summary.DataTableWithClasses <- function(self) {
#  message("\t\t\n giving the summary of the created object")
  cat("DataTableWithClasses\n")
  cat("\tData type         \t", self$dataType, "\n")
  # cat("\tVariables type         \t", x$variablesType, "\n")
  cat("\tNumber of genes   \t", self$nbGenes, "\n")
  cat("\tNumber of samples \t", self$nbSamples, "\n")
  cat("\tNumber of classes \t", self$nbClasses, "\n")
  cat("\tClass properties\n")
  print(self$classProperties)
  cat("\n")
  NextMethod("summary", self)
}


#' @title finish summary
#' @author Mustafa AubElQumsan and Jacques van Helden
#' @export
summary.default <- function(self) {
   cat("")
}

#' @export
print.DataTableWithClasses <- function(x) {
  summary(x)
}




#' @title build the attributes of an object of class DataTableWithClasses based on the count table and pheno table.
#' @author Jacques van Helden and Mustafa AbuElQumsan
#' @description Derive the attributes of an object from the count table and pheno table.
#' This includes the definition of sample and gene names, as well as the association of
#' classes to samples, based on the classColumn(s) specified in the parameters.
#' @parameters self an object of class DataTableWithClasses
#' @return the same object with added fields to describe sample classes and derived
#' attributes. Sample classes are extracted from the phenoTable, using column(s)
#' specified in self$classColumn.
#'
#' @export
buildAttributes.DataTableWithClasses <- function(self) {
  message("\tBuilding class-specific attributes for DataTableWithClasses\t", self$ID)

  ## Check rows of pheno table
  if (nrow(self$phenoTable) != ncol(self$dataTable)) {
    stop("DataTableWithClasses(): inconsistent dimensions of phenoTable (",
         nrow(phenoTable), " rows) and dataTable (",
         ncol(dataTable), " columns).")
  }

  ## Check if the pheno table contains a column corresponding to the indication
  if (sum(!self$classColumn %in% names(self$phenoTable)) > 1) {
    stop("\n\tMissing column(s) in the pheno table ", paste(collapse=", ", setdiff(classColumn, names(phenoTable))),
         "\n\tColumns found in the pheno table: ", paste(collapse=", ", names(self$phenoTable)))
  }

  ## Count table-derived paraemters
  self$nbSamples = ncol(self$dataTable)
  self$nbGenes = nrow(self$dataTable)
  self$sampleNames = colnames(self$dataTable)
  self$geneNames = rownames(self$dataTable)

  ################################################################
  ## Specify sample classes (classLabels) and all the derived attributes
  ## (classNames, nbClasses, ...) by extracting information about
  ## specified class columns
  if (is.null(self$classColumn) || (length(self$classColumn) < 1)) {
    stop("classColumn must be defined. ")
  } else if (length(self$classColumn) == 1) {
    self$classLabels <-  as.vector(self$phenoTable[, self$classColumn])
  } else {
    ## Combine several columns to establish the classLabels
    self$classLabels <- apply(self$phenoTable[, self$classColumn], 1, paste, collapse="_")
  }
  # table(classLabels)
  self$classNames <- unique(sort(self$classLabels))
  self$nbClasses <- length(self$classNames)

  ## Build a table with class properties (size, color, ...)
  self$classProperties <- (as.data.frame.table(table(self$classLabels)))
  colnames(self$classProperties) <- c("Class", "nbSamples")
  rownames(self$classProperties) <- as.vector(self$classProperties$Class)

  ## Build a vector with the number of samples per class
  self$samplesPerClass <- as.vector(as.matrix(self$classProperties$nbSamples))
  names(self$samplesPerClass) <- as.vector(as.matrix(self$classProperties$Class))

  ## Relative frequencies of individuals per class
  self$classFrequencies <- self$samplesPerClass / sum(self$samplesPerClass)
  self$randExpectedHitRate = self$classFrequencies %*% self$classFrequencies
  self$randExpectedMisclassificationRate = 1 - self$randExpectedHitRate

  ## Define class colors
  if (is.null(self$classColors)) {
    self$classColors <- rainbow(n  = self$nbClasses)
    names(self$classColors) <- self$classNames
    #classColors <- unlist(classColors)
  } else {
    if (class(self$classColors) == "list") {
      ## Convert list (like the one parsed  from the YAML file) to a named vector
      self$classColors <- unlist(self$classColors)
    }

    ## Assign automatic colors (numbers starting from 1) for classes with no defined color in the parameters
    checkedColors <- self$classColors[self$classNames]
    names(checkedColors) <- self$classNames
    missing.color <- is.na(checkedColors)
    if (sum(missing.color) > 0) {
      ## Pick up missing color in a rainbow with n colors = nbClasses
      checkedColors[missing.color] <- rainbow(n = self$nbClasses)[missing.color]
      self$classColors <- checkedColors
    }


    ##    self$classNames %in% names(self$classColors)
  }

  ## Assign class-specific colors if they were defined in the parameters
  self$classProperties$color <- self$classColors[rownames(self$classProperties)]


  ## Assign colors to samples
  self$sampleColors <- self$classColors[self$classLabels]
  names(self$sampleColors) <- self$sampleNames

  #### Check consistency of  parameters ####


  ## Check sample names
  if (length(self$sampleNames) != ncol(self$dataTable)) {
    stop("DataTableWithClasses(): inconsistent dimensions of sampleNames (",
         length(self$sampleNames), " names) and dataTable (",
         ncol(self$dataTable), " columns).")
  }

  ## Check gene names
  if (length(self$geneNames) != nrow(self$dataTable)) {
    stop("DataTableWithClasses(): inconsistent dimensions of geneNames (",
         length(self$geneNames), " names) and dataTable (",
         nrow(self$dataTable), " columns).")
  }

  self <- NextMethod("buildAttributes", self)
  #  message("trainIndices length : ", length(self$trainTestProperties$trainIndices))
  #  message("\t\treturning from buildAttributes.DataTableWithClasses()")
  return (self)
}



