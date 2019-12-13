require("RNAseqMVA")

##### Path of the YAML-formatted configuration file ####
configFile <- "misc/00_project_parameters.yml"

#### Read parameters from the yaml-formatted config file ####
## These parameters will then be used to overwrite the default parameters.
message.with.time("Loading parameters from YAM file ", configFile)
project.parameters <- yaml.load_file(configFile)


#### Read arguments passed from the command line ####
args = commandArgs(trailingOnly=TRUE)
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

#### RecountIDs with problems (TO INVESTIGATE / DEBUG LATER) ####
recountIDs.with.problems <- c("SRP003611" = "the number of samples drops after run merging",
                              "SRP039694" = '	Building class-specific attributes for DataTableWithClasses	SRP039694
Error in `colnames<-`(`*tmp*`, value = c("Class", "nbSamples")) :
  "names" attribute [2] must be the same length as the vector [1]',
                              "SRP008976" = "",
                              "SRP006574" = "	Building pheno table\nError in rbind(deparse.level, ...) : numbers of columns of arguments do not match",
                              "SRP042161" = "class column not specified",
                              "SRP041736" = "We cannot analyse this dataset because the pheno table does not contain any info about the sample classes
class column not specified (NA)
Instead of complaning because of this, the program crashes with the following message:

TO DO: detect such cases and stop with explicit message before any analysis.

2018-04-07_093304		Creating object of class DataTableWithClasses
	Building attributes for object of class DataTableWithClasses
                              Building class-specific attributes for DataTableWithClasses	SRP041736
                              Error: subscript contains invalid names
                              ")

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

## Note Jv (2019-10-29) I moved initParallelComputing() to a separate script
##01C_init_parallel_computing.R
## because we sometimes want to reload the parameters without re-creating new
## slots for parallel computing

#### END OF SCRIPT ####

