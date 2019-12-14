#' @title Intialise parameters and directories for a give recount ID.
#' @author Jacques van Helden and Mustafa AbuElQumsan
#' @description This function ensures the consistency of the parameters for a giveb RNA-seq experiment specified by its recount ID.
#' (1) Merge the global parameters with experiment-specific parameters
#' (description, short label, ...).
#' (2) Check that some required parameters are well defined.
#' (3) Check the output directory (if not specified in config, define them with default paths) and create them if required.
#' @param recountID must be a valid ID row Recount2 database (SRP...)
#' @param project.parameters project parameters (read from the YAML-formatted file)
#' @return all project.parameters for the related recount IDn as specified in the YAML configuration file.
#' @export
initRecountID <- function(recountID, project.parameters) {
  message.with.time("Loading parameters from YAM file ", configFile)


  ## Load default parameters for each new recountID
  ## (was previously parsed from the YAML file)
  parameters <- c(project.parameters$default, project.parameters$global)


  ## Specify the current recountID in parameters
  parameters$recountID <- recountID

  ## Get the feature type (gene or transcript)
  if (is.null(parameters$feature)) {
    parameters$feature <- "gene"
    message("Beware: YAML configuration file did not specify the feature type -> set to gene")
  }
  message("\tFeature type\t", parameters$feature)

  ## studyPath, i.e. path where the data downloaded from Recount2 is stored
  ## Note: we store the different feature types (gene, transcript)
  ## in the same folder, to mirror  Recount2 server
  parameters$studyPath <- file.path(parameters$dir$workspace, "data", recountID)



  ## Overwrite default parameters wih project-specific parameters
  specific.parameters <- project.parameters[[recountID]]
  if (is.null(specific.parameters)) {
    message("\tNo specific parameters for recount ID ", recountID)
    message("\tUsing generic parameters from the yaml file. ")
  } else {
    message("\tUsing specific parameters specfied in yaml file for recount ID ", recountID)
    parameters[names(specific.parameters)] <- project.parameters[[recountID]]
    names(parameters$data.types) <- parameters$data.types
    names(parameters$variables.type) <- parameters$variables.type
  }


  #### Check some required parameters ####


  ## Classifiers must be specified
  if (is.null(parameters$classifiers)) {
    stop("The configuration file must contain a parameter named 'classifiers'.")
  }

  ## The list of classifiers cannot be empty
  if (length(parameters$classifiers) < 1) {
    stop("The configuration file must specify at least one classifier.")
  }

  ## Convert list-formatted class colors to named vector (YAML does not allow to specify named vectors)
  if (!is.null(parameters$classColors)) {
    if (class(parameters$classColors) == "list") {
      parameters$classColors <- unlist(parameters$classColors)
    }
  }

  ##### Output directories ####
  if (is.null(parameters$dir)) {
    stop("Output directories (main, workspace) must be specified in the config file. ")
  }


  ## Check required directories
  if (is.null(parameters$dir$main)) {
    stop("Main directory must be specified in the config file. ")
  }

  ## Workspace
  if (is.null(parameters$dir$workspace)) {
    stop("Workspace directory must be specified in the config file. ")
  }

  ## Reset the result dir for this particular study

  if (is.null(parameters$dir$results)) {
    parameters$dir$results <- file.path(
      parameters$dir$workspace,
      "results",
      paste0(recountID, "_",
             parameters$feature))
  } else {
    parameters$dir$results <- file.path(
      parameters$dir$results,
      paste0(recountID, "_", parameters$feature))
  }
  dir.create(parameters$dir$results, showWarnings = FALSE, recursive = TRUE)
  message("\t\tresults\t", parameters$dir$results)

  ## Sub-directory to save the generic figures of this study case (for example, comparisons between study cases)
  if (is.null(parameters$dir$figures)) {
    parameters$dir$figures <- file.path(
      parameters$dir$results,
      "figures")
  }
  dir.create(path = parameters$dir$figures, recursive = TRUE, showWarnings = FALSE)
  message("\t\tfigures\t", parameters$dir$figures)

  ## Sub-directory to save the tab-separated value (TSV) files
  if (is.null(parameters$dir$tsv)) {
    parameters$dir$tsv <- file.path(
      parameters$dir$results,
      "TSV")
  }
  dir.create(path = parameters$dir$tsv, recursive = TRUE, showWarnings = FALSE)
  message("\t\tTSV\t", parameters$dir$tsv)


  #### Define the classifier-specific directories where tables and figures will be stored. ####
  ## One sub-directory per classifer, with separate subdirectories for tables and figures.
  if (is.null(parameters$dir$classifiers)) {
    parameters$dir$classifiers <- file.path(
      parameters$dir$results,
      parameters$classifiers)
  }
  names(parameters$dir$classifiers) <- parameters$classifiers
  message("\tclassifiers\n\t\t", paste(collapse="\n\t\t", parameters$dir$classifiers))

  if (is.null(parameters$dir$classifier_tables)) {
    parameters$dir$classifier_tables <- file.path(
      parameters$dir$classifiers, "tables")
  }
  names(parameters$dir$classifier_tables) <- parameters$classifiers
  message("\tclassifier_tables\n\t\t", paste(collapse = "\n\t\t", parameters$dir$classifier_tables))

  if (is.null(parameters$dir$classifier_figures)) {
    parameters$dir$classifier_figures <- file.path(
      parameters$dir$classifiers, "figures")
  }
  names(parameters$dir$classifier_figures) <- parameters$classifiers
  message("\tclassifier_figures\n\t\t", paste(collapse = "\n\t\t", parameters$dir$classifier_figures))

  ## Create all the recountID-specific sub-directories
  for (dir in c(parameters$dir$classifiers,
                parameters$dir$classifier_tables,
                parameters$dir$classifier_figures
                )) {
    dir.create(dir, showWarnings = F, recursive = T)
  } # end loop over the dir

  #### TO CHECK LATER: DO wE STILL NEED THESE VARIABLES ??? ####

  ## Directory for impact of Normalization and log2 into counts (and the study of its impact)
  if (is.null(parameters$dir$NormalizationImpact)) {
    parameters$dir$NormalizationImpact  <- file.path(
      parameters$dir$results,
      paste("impact_of_normalisation", sep = ""))
  }
  dir.create(parameters$dir$NormalizationImpact, showWarnings = F, recursive = T)
  message("\t\tNormalizationImpact\t", parameters$dir$NormalizationImpact)

  ## Directory for the visualization of Principal component for counts (and the study of its impact)
  if (is.null(parameters$dir$PCviz)) {
    parameters$dir$PCviz <- file.path(
      parameters$dir$results,
      paste( "visualization_of_PCs", sep = ""))
  }
  dir.create(parameters$dir$PCviz, showWarnings = F, recursive = T)
  message("\t\tPCviz\t", parameters$dir$PCviz)

  # View(parameters)
  return(parameters)
}
