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

  if (is.null(parameters$dir$main)) {
    stop("Main directory must be specified in the config file. ")
  }

  # ## Directory with R scripts
  # if (is.null(parameters$dir$scripts)) {
  #   parameters$dir$scripts <- file.path(dir.main, "R")
  # }

  return(parameters)
}
