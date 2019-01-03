#### Reload parameters and update their values in each study case  ####

## If requested, reset the parameters for all the study cases
## This is used to re-run the analyses on each study case after
## having changed some parameters in the yaml-specific configuration file
if (!is.null(project.parameters$global$reload.parameters)
    && project.parameters$global$reload.parameters) {
  message.with.time("Reloading parameters from config file\t", configFile)
  project.parameters <- yaml.load_file(configFile)
  project.parameters <- initParallelComputing(project.parameters)
  if (exists("studyCases")) {
    for (recountID in names(studyCases)) {
      parameters <- initRecountID(recountID, project.parameters)
      studyCases[[recountID]]$parameters <- parameters
      for (datasetName in names(studyCases[[recountID]]$datasetsForTest)) {
        studyCases[[recountID]]$datasetsForTest[[datasetName]]$parameters <- parameters
      }
    }
  }
}

message.with.time("Finished running script\t", "misc/01c_reload_parameters.R")

