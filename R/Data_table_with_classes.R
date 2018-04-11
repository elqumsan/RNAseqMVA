
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

