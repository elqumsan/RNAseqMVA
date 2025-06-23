#### Load parameters from a YAML-formatted configuration file ####
#### This script is used to load the parameters for the RNAseqMVA project
#### It is used in the main script (RNAseqMVA.R) and in the scripts for each recount ID

require("RNAseqMVA")

##### Path of the YAML-formatted configuration file ####
configFile <- "misc/00_project_parameters.yml"

#### Read parameters from the yaml-formatted config file ####
## These parameters will then be used to overwrite the default parameters.
message.with.time("Loading parameters from YAM file ", configFile)
project.parameters <- yaml.load_file(configFile)

#### Read arguments passed from the command line ####
args = commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  project.parameters$global[["selected_recount_ids"]] <- args[1]
  message("Overwriting selected_recount_ids parameter with command-line argument ", args[1])
}
if (length(args) >= 2) {
  message("Overwriting config feature parameter with command-line argument ", args[2])
  project.parameters$global[["feature"]] <- args[2]
}


#### Specify selected recount IDs ####
selectedRecountIDs <- project.parameters$global[["selected_recount_ids"]]


## Get the list of recountIDs with parameters specified in the YAML file
recountIDs <- grep(pattern = "^SRP", x = names(project.parameters), value = TRUE)
message("\tYAML config file contains ", length(recountIDs)," recount IDs.")
message("\tSelected ", length(selectedRecountIDs)," recount IDs: ", paste(collapse = "; ", selectedRecountIDs))

## Feature type
featureType <- project.parameters$global$feature
message("\tFeature type: ", featureType)


#### Initialise global directories ####

## Directory for cross-study case comparisons
if (is.null(project.parameters$global$dir$figures)) {
  project.parameters$global$dir$figures <- file.path(project.parameters$global$dir$results, "figures")
}

for (d in names(project.parameters$global$dir)) {
  dir <- project.parameters$global$dir[[d]]
  if (!dir.exists(dir)) {
    message("Creating ", d, " directory\t", dir)
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  }
}

# View(project.parameters)

#### END OF SCRIPT ####

