#' @title Iintialise parameters and directories for a give recount ID.
#' @author Jacques van Helden and Mustafa AbuElQumsan
#' @param recountID must be a valid ID rom Recount2 database (SRP...)
#' @param configFile YAML-formatted configuration file specifying global and recountID-specific parameters
#' @export
initRecountID <- function(recountID, configFile) {
  message.with.time("Loading parameters from YAM file ", configFile)

  project.parameters <- yaml.load_file(configFile)

  ## Load default parameters for each new recountID
  ## (was previously parsed from the YAML file)
  parameters <- project.parameters$default

  ## Specify the current recountID in parameters
  parameters$recountID <- recountID

  ## Overwrite default parameters wih project-specific parameters
  selected.parameters <- project.parameters[[recountID]]
  if (is.null(selected.parameters)) {
    message("No specific parameters for recount ID ", recountID)
    message("Using generic parameters from the yaml file. ")
  } else {
    message("Using specific parameters specfied in yaml file for recount ID ", recountID)
    parameters[names(selected.parameters)] <- project.parameters[[recountID]]
    names(parameters$data.types)<-parameters$data.types
    names(parameters$variables.type)<-parameters$variables.type
  }

  ## Convert list-formatted class colors to named vector (YAML does not allow to specify named vectors)
  if (!is.null(parameters$classColors)) {
    if (class(parameters$classColors) == "list") {
      parameters$classColors <- unlist(parameters$classColors)
    }
  }

  ## Specify directories
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


  parameters$dir$results <- file.path(parameters$dir$workspace, "results", recountID)


  # ## Directory with R scripts
  # if (is.null(parameters$dir$scripts)) {
  #   parameters$dir$scripts <- file.path(dir.main, "R")
  # }


  ## Result directory
  parameters$dir$results <- file.path(parameters$dir$workspace, "results", parameters$recountID)

  ## Sub-directory to save the tab-separated value (TSV) files
  parameters$dir$tsv <- file.path(parameters$dir$results, "TSV")
  dir.create(path = tsv.dir, recursive = TRUE, showWarnings = FALSE)


  return(parameters)
}
