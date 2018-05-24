#' @title Intialise parameters and directories for a give recount ID.
#' @author Jacques van Helden and Mustafa AbuElQumsan
#' @param recountID must be a valid ID row Recount2 database (SRP...)
#' @param project.parameters project parameters (read from the YAML-formatted file)
#' @return all project.parameters for the related recount ID.
#' \itemize {
#'      \item global: ## Global parameters are those who are the samee across all study cases
#'      \item  dir: ## Directories
#'      \item   main: "~/RNAseqMVA"
#'      \item  workspace: "~/RNAseqMVA_workspace"
#'      \item results: "~/RNAseqMVA_workspace/results"
#'      \item TSV: "~/RNAseqMVA_workspace/data"
#'      \item memoryImages: "~/RNAseqMVA_workspace/memory_images"

#'      \item perm.prefix: "permLabels"
#'      \item jobs: "none" ## Maximal number of jobs to run in parallel. Supported: integers > 0, "auto", or "none"
#  iterations: 50, ## Number of iterations for the classiifers
#'       \item iterations: 10 ## Number of iterations for the classiifers

#'       \item reload: FALSE
#'       \item compute: TRUE ## If FALSE, do not run the heavy computations, just generate the pictures and save tables
#'       \item save.tables: TRUE ## If TRUE, save tab-delimited files with the results
#'       \item save.image: FALSE
#' draw.plot : TRUE
#'        \item classifiers: "svm"
# classifiers: ["svm", "knn","rf"]
# classifiers: ["knn","rf","svm"] # TO ADD LATER: "lda", "qda" + the one suggested by Aitor

#'       \item data.types: ["countsPerRun", "originalCounts", "filtered", "norm", "log2norm", "log2normPCs", "DEG", "V.importance"]
#'       \item data.types.to.test: ["log2normPCs"]
# data.types.to.test: ["filtered", "norm", "log2norm", "log2normPCs"]

# data.types.to.test: ["log2norm"]
#'       \item variables.type: ["all", "top_ordered"]
#'       \item deg.methods: ["DESeq2", "edgeR"]
## Note: for the number of variables I choose a regular spacing around the total number of samples in order to detect overfitting effects
#'       \item nb.variables: [3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 300, 400, 500, 1000, 2000, 5000, 10000] ## Number of variables for variable ordering test
#'       \item trainingProportion: [0.66] # ratio of spliting the data set to training and testing sets
#'       \item identicalTrainTest: TRUE
#'       \item permute: [FALSE, TRUE]
#'       \item verbose: TRUE

#'       \item knn:
#'       \item   k: 10
#'       \item k_values: [3, 5, 7, 10, 15]

#'       \item svm:
#'       \item kernel: "linear"
#'       \item kernel_values: ["linear", "polynomial", "sigmoid"] ## NOT WORKING YET, TO BE TESTED LATER: "radial basis"
#'       \item scale: FALSE
#'       \item type: "C-classification"
# type: ["C-classification", "nu-classification", "one-classification"]
#'       \item gamma: 1
#'       \item cost: 100

## Default values for study case-specific parameters
## (can then be over-written by recountID-specific values as shown below)
#'      \item default:
#  classColumn: "tissue" ## Characteristics to be used as class label
#'      \item mergeRuns: TRUE ## Whether or not to merge runs per sample
#'      \item sampleIdColumn: "geo_accession"

## Filtering options
#'      \item minSamplesPerClass: 15 ## THERE IS A BUG TO FIX: if we set this parameter to 10 (for example) the class labels must be adapted by suppressing the
#'      \item nearZeroVarFilter: TRUE
#'      \item na.rm: TRUE

#' }

#' @export
initRecountID <- function(recountID, project.parameters) {
  message.with.time("Loading parameters from YAM file ", configFile)


  ## Load default parameters for each new recountID
  ## (was previously parsed from the YAML file)
  parameters <- c(project.parameters$default, project.parameters$global)

  ## Specify the current recountID in parameters
  parameters$recountID <- recountID
  parameters$studyPath <- file.path(parameters$dir$workspace, "data", recountID)


  ## Overwrite default parameters wih project-specific parameters
  selected.parameters <- project.parameters[[recountID]]
  if (is.null(selected.parameters)) {
    message("No specific parameters for recount ID ", recountID)
    message("Using generic parameters from the yaml file. ")
  } else {
    message("Using specific parameters specfied in yaml file for recount ID ", recountID)
    parameters[names(selected.parameters)] <- project.parameters[[recountID]]
    names(parameters$data.types) <- parameters$data.types
    names(parameters$variables.type) <- parameters$variables.type
  }


  #### Check some required parameters ####

  ## Classifiers must be specified
  if (is.null(parameters$classifiers)) {
    stop("The configuration file must contain a parameter named 'classifiers'.")
  }

  if (length(parameters$classifiers) < 1) {
    stop("The configuration file must specify at least one classifier.")
  }

  ## Convert list-formatted class colors to named vector (YAML does not allow to specify named vectors)
  if (!is.null(parameters$classColors)) {
    if (class(parameters$classColors) == "list") {
      parameters$classColors <- unlist(parameters$classColors)
    }
  }

  ##### Specify directories ####
  if (is.null(parameters$dir)) {
    parameters$dir <- list()
  }

  ## Check required directories
  if (is.null(parameters$dir$main)) {
    stop("Main directory must be specified in the config file. ")
  }

  ## Workspace
  if (is.null(parameters$dir$workspace)) {
    stop("Workspace directory must be specified in the config file. ")
  }

  ## Result directory
  parameters$dir$results <- file.path(parameters$dir$workspace, "results")
  dir.create(parameters$dir$results, showWarnings = FALSE, recursive = TRUE)

  # ## Directory with R scripts
  # if (is.null(parameters$dir$scripts)) {
  #   parameters$dir$scripts <- file.path(dir.main, "R")
  # }


  ## Result directory
  #parameters$dir$results <- file.path(parameters$dir$workspace, "results", parameters$recountID)

  ## Sub-directory to save the tab-separated value (TSV) files
  parameters$dir$tsv <- file.path(parameters$dir$results,recountID, "TSV")
  dir.create(path = parameters$dir$tsv, recursive = TRUE, showWarnings = FALSE)

  ## Define the classifier-specific directories where tables and figures will be stored.
  ## One sub-directory per classifer, with separate subdirectories for tables and figures.

  parameters$dir$classifiers <- file.path( parameters$dir$results,recountID,parameters$classifiers)
  names(parameters$dir$classifiers) <- parameters$classifiers

  parameters$dir$tables <- file.path(parameters$dir$classifiers, "tables")
  names(parameters$dir$tables) <- parameters$classifiers

  parameters$dir$figures <- file.path(parameters$dir$classifiers, "figures")
  names(parameters$dir$figures) <- parameters$classifiers

  parameters$dir$figuresDetail <- file.path(parameters$dir$figures, "detailFigures")
  names(parameters$dir$figuresDetail) <- parameters$classifiers

  parameters$dir$tablesDetail <- file.path(parameters$dir$tables, "detailTables")
  names(parameters$dir$tablesDetail) <- parameters$classifiers

  ## Create all the recountID-specific sub-directories
  for (dir in c(parameters$dir$classifiers, parameters$dir$tables, parameters$dir$figures, parameters$dir$figuresDetail, parameters$dir$tablesDetail)) {
    dir.create(dir, showWarnings = F, recursive = T)
  } # end loop over the dir

  ################################################################
  ## TO CHECK LATER: DO wE STILL NEED THESE VARIABLES ???

  ## Directory for impact of Normalization and log2 into counts (and the study of its impact)
  parameters$dir$NormalizationImpact <- file.path(parameters$dir$results, recountID , paste("impact_of_normalisation_and_log2", sep = ""))
  dir.create(parameters$dir$NormalizationImpact, showWarnings = F, recursive = T)

  ## Directory for the visualization of Principal component for counts (and the study of its impact)
  parameters$dir$PCviz <- file.path(parameters$dir$results, recountID , paste( "visualization_of_PCs", sep = ""))
  dir.create(parameters$dir$PCviz, showWarnings = F, recursive = T)

  # View(parameters)
  return(parameters)
}
