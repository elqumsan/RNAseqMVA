#### Load parameters from a YAML-formatted configuration file ####
#### This script is used to load the parameters for the RNAseqMVA project
#### It is used in the main script (RNAseqMVA.R) and in the scripts for each recount ID

require("RNAseqMVA")
require("optparse")

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

## Read the command-line arguments passed to Rscript with the optparse package
## This is used to overwrite the parameters specified in the YAML file
## This is useful for running the script with different parameters without modifying the YAML file
## For example, to run the script with a different Recount ID and feature type
##    Rscript RNAseqMVA.R --recountID SRP048759 --feature gene

# Define options
option_list = list(
  make_option(c("-r", "--recountID"), type="character", default=NULL,
              help="Recount ID (e.g. SRP056295)", metavar="character"),
  make_option(c("-f", "--feature"), type="character", default="gene",
              help="Feature type: gene or transcript [default: %default]", metavar="character")
)

# Parse options
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

## Update recountID if specified on the command-line arguments

if (!is.null(opt$recountID)) {

  # Example use
  message("Recount ID specified in the command-line arguments: ", opt$recountID, "\n")


  ## Check that the yaml config file contains metadata for the command-line specified recountID
  ## If not, stop the script with an error message

  if (is.null(project.parameters[[opt$recountID]])) {
    stop("Recount ID '", opt$recountID, "' is not specified in the YAML configuration file. ",
         "Please check the config file or specify a valid recount ID.")
  }

  selectedRecountIDs[1] <- opt$recountID  # Update the selected recount IDs with the command-line specified recount ID
}

## Update feature type if specified on the command-line arguments

if (!is.null(opt$feature)) {

  # Validate feature argument
  valid_features <- c("gene", "transcript")
  if (!(opt$feature %in% valid_features)) {
    stop("Invalid feature type '", opt$feature,
         "'. Valid values are: ", paste(valid_features, collapse = ", "))
  }

  # Handle the feature type
  message("\tFeature type specified in the command-line arguments: ", opt$feature, "\n")
  project.parameters$global$feature <- opt$feature
}


#### END OF SCRIPT ####

