#' @title Intialise parameters and directories for a give recount ID.
#' @author Jacques van Helden and Mustafa AbuElQumsan
#' @param recountID must be a valid ID row Recount2 database (SRP...)
#' @param project.parameters project parameters (read from the YAML-formatted file)
#' @export
initRecountID <- function(recountID, project.parameters) {
  message.with.time("Loading parameters from YAM file ", configFile)


  ## Load default parameters for each new recountID
  ## (was previously parsed from the YAML file)
  parameters <- project.parameters$default

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
