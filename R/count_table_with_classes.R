#' @title constructor of the countTableWithClasses class
#' @author Mustafa AbuElQumsan and Jacques van helden
#' @description This class contains count tables
#' (typically used to represent expression profiles in RNA-seq analysis)
#' documented with a Pheno table (description of each sample) and completed
#' with some additional fields for convenience: sample names, gene names, ...
#' @param countTable a data.frame with one row per feature (e.g. gene) and
#' one column per sample.
#' @param phenoTable a data frame describing each sample: one row per sample and one column per attribute of a sample.
#' @param classColumn a vector indicating one or several columns of the pheno table that will define the sample class labels. If umtiple columns are specified, they will be concatenated to build sample labels
#' @param ID=parameters$recountID identifier associated to the count table (by default, the RecountID, but can be specified with custom identifiers)
#' @param sampleNames=colnames(countTable) sample names (by default, automatically taken from the column names of the count table)
#' @param geneNames=rownames(countTable) gene names (by default, automatically taken from the row names of the count table)
#' @param dataType="raw counts" data type, free text (e.g. raw counts, log2-transformed counts, log2 normalised counts  ...).
#'
#' @export

countTableWithClasses <- function(countTable,
                                  phenoTable,
                                  classColumn,
                                  classColors = parameters$classColor,
                                  ID = parameters$recountID,
                                  sampleNames = colnames(countTable),
                                  geneNames = rownames(countTable),
                                  dataType = "raw counts"
) {

  ## Built a list from the input parameters
message.with.time("\tCreate object has the countTableWithClasses attribute" )

  object <- structure(
    list(
      ID = ID,
      countTable = countTable,
      phenoTable = phenoTable,
      sampleNames = sampleNames,
      nbSamples = ncol(countTable),
      geneNames = geneNames,
      nbGenes = nrow(countTable),
      dataType = dataType
    ),
    class="countTableWithClasses")
  # names(object)
  # class(object)
  # attributes(object)


  #### Check consistency of input parameters ####

  ## Check rows of pheno table
  if (nrow(object$phenoTable) != ncol(object$countTable)) {
    stop("countTableWithClasses(): inconsistent dimensions of phenoTable (",
         nrow(phenoTable), " rows) and countTable (",
         ncol(countTable), " columns).")
  }

  ## Check sample names
  if (length(object$sampleNames) != ncol(countTable)) {
    stop("countTableWithClasses(): inconsistent dimensions of sampleNames (",
         length(object$sampleNames), " names) and countTable (",
         ncol(countTable), " columns).")
  }

  ## Check gene names
  if (length(object$geneNames) != nrow(countTable)) {
    stop("countTableWithClasses(): inconsistent dimensions of geneNames (",
         length(object$geneNames), " names) and countTable (",
         nrow(countTable), " columns).")
  }


  ################################################################
  ## Specify sample classes (classLabels) by extracting information about specified class columns
  if (is.null(classColumn) || (length(classColumn) < 1)) {
    stop("classColumn must be defined. ")
  } else if (length(classColumn) == 1) {
    object$classLabels <-  as.vector(object$phenoTable[, classColumn])
  } else {
    ## Combine several columns to establish the classLabels
    object$classLabels <- apply(object$phenoTable[, classColumn], 1, paste, collapse="_")
  }
  # table(classLabels)
  object$classNames <- unique(sort(object$classLabels))
  object$nbClasses <- length(object$classNames)

  ## Build a table with class properties (size, color, ...)
  object$classProperties <- (as.data.frame.table(table(object$classLabels)))
  colnames(object$classProperties) <- c("Class", "nbSamples")

  ## Build a vector with the number of samples per class
  object$samplesPerClass <- as.vector(as.matrix(object$classProperties$nbSamples))
  names(object$samplesPerClass) <- as.vector(as.matrix(object$classProperties$Class))

  ## Define class colors
  classColors <-1:length(object$classNames)
  names(classColors)<- object$classNames
  #classColors <- unlist(classColors)
  object$classProperties$color<- classColors

  ## Assign colors to samples
  object$sampleColors <- classColors[object$classLabels]
  names(object$sampleColors) <- object$sampleNames

  # ## Class colors may be defined in the yaml parameters
  # if (!is.null(parameters$classColor)) {
  #   #  parameters$classColor[["astrocytes"]]
  #   ## Convert the yaml-imported list into a named vector
    # convert list to a vector

    #classColors <- 1:length(loaded$filtered$classColors)
    #names(loaded$filtered$classColors) <- loaded$filtered$classNames


    #names(sampleColors) <- rownames(loaded$filtered$phenoTable)
    # names(classColors) # check vector names

    # classColors["astrocytes"]
  #
  # } else {
  #   # # loaded$filtered$classNames <- unique(loaded$filtered$classLabels)
  #   #  classColors <- 1:length(classNames)
  #   #  names(classColors) <- classNames
  #   stop("Don't haveing any class Color, you should go back yml file to revise it...")
  #
  # }


  message("\t\tfinishing from creating the object with countTablewithClasses attribute")
  return(object)
}


summary.countTableWithClasses <- function(x) {
#  message("\t\t\n giving the summary of the created object")
  cat("countTableWithClasses\n")
  cat("\tData type         \t", x$dataType, "\n")
  cat("\tNumber of genes   \t", x$nbGenes, "\n")
  cat("\tNumber of samples \t", x$nbSamples, "\n")
  cat("\tNumber of classes \t", x$nbClasses, "\n")
  cat("\tClass properties\n")
  print(x$classProperties)
  cat("\n")
}

print.countTableWithClasses <- function(x) {
  summary.countTableWithClasses(x)
}

## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ##
## TEMPORARILY HERE
## We define a method for stratified selection of training and testing sets,
## In a second time we will creat a new class names TrainTestCounttable, and this
## method will be attached to this class.

#' @title Random sampling of n training sets in a CountTableWithClasses.
#' @author Jacques van Helden and Mustafa AbuElQumsan
#' @description Select n random subsets for training. among the biological samples from a CountTableWithClasses.
#' @param self an object of the class CountTableWithClasses
#' @paraam stratified=TRUE if true, sampling is done in each class separately in order to preserve the relative frequencies of classes in training and testing sets.
#' @param iterations=parameters$iterations number of  train/test  iterations, which defines the number of independent sampled subsets
#' @param trainingProportion=parameters$trainingProportion proportion of samples to sample for each training set
#'
#' @export

selectTrainingSets <- function(countTable,
                               stratified=TRUE,
                               iterations = parameters$iterations,
                               trainingProportion = parameters$trainingProportion) {
  message.with.time("Selecting ", iterations, " training sets, with training proportion = ", trainingProportion)

  #### Check validity of the paraemters ####

  ## Check the class of input object
  if (!is(countTable, "countTableWithClasses")) {
    stop("selectStratifiedTrainingSets(): countTable parameter should belong to class countTableWithClasses. ")
  }
  ##  STRANGE: THIS RETURNS FALSE WHEREAS IT SHOULD B TRUE
  # isClass("countTableWithClasses")


  ## Trainng Proportion
  if ((trainingProportion < 0) || (trainingProportion > 1)) {
    stop("Training proportion must be a real number comprised between 0 and 1")
  }

  ## Instantiate the list with training indices
  trainIndices <- list()

  if (stratified) {
    ## Get class sizes
    trainSizePerClass <- round(countTable$samplesPerClass * trainingProportion)
    message("Stratified sampling among classes")
    print(as.data.frame(trainSizePerClass))
    for (i in 1:parameters$iterations) {
      trainIndices[[i]] <- vector()
      # c <- 2
      for (c in 1:countTable$nbClasses) {
        currentClass <- countTable$classNames[[c]]
        classSamples <- which (countTable$classLabels == currentClass)
        classTrain <- sample(classSamples, size = trainSizePerClass[[currentClass]], replace = FALSE)
        trainIndices[[i]] <- append(trainIndices[[i]], classTrain)
        ## Check that the stratification  was correct
        ## table(countTable$classLabels[trainIndices[[i]]]) == trainSizePerClass
      }
    }
  } else {
    ## Sample the training sets irrespective of class membership
    n <- countTable$nbSamples
    train.size <- round(trainingProportion * n)
    message("Class-independent sampling of training sets")
    for (i in 1:parameters$iterations) {
      trainIndices [[i]] <- sample(1:n, size = train.size, replace = FALSE)
    #  View(as.data.frame.list(trainIndices))
    }
  }


  ## Select testIndices as the complement of train indices
  testIndices <- list()
  for (i in 1:parameters$iterations) {
    ## MUTSAFA: DO IT
  }


  ## Add the attributes
  countTable$iterations <- iterations
  countTable$trainingProportion <- trainingProportion
  countTable$trainIndices <- trainIndices



  # return(countTable)
  message.with.time("Stratified training set selection done")
}
