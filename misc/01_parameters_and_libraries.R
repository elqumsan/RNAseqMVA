require("RNAseqMVA")

message.with.time("Loading required libraries")

require("yaml")
require("roxygen2") ## MUSTAFA: do we really required roxygen2 and devtools to run the scripts ? They are required for building the package, but maybe not for running the analyses
require("devtools") ## MUSTAFA: do we really required roxygen2 and devtools to run the scripts ? They are required for building the package, but maybe not for running the analyses

# loading required libraries
requiredCRAN <- c('devtools', 'class', "randomForest","broom", "roxygen2","scatterplot3d","e1071")
RequiredCRANPackages(requiredCRAN)

## JvH: Mustafa, please add the other required packages, in particular recount
requiredBioconductor <- c("recount")
RequiredBioconductorPackages(requiredBioconductor)

yaml.file <- "~/RNAseqMVA/misc/00_project_parameters.yml"
message.with.time("Loading parameters from YAM file ", yaml.file)

## Read project-specific parameters from a yaml-formatted file.
## These parameters will then be used to overwrite the default parameters.
project.parameters <- yaml.load_file(yaml.file)


## Get all recount IDs
recountIDs <- grep(pattern = "^SRP", x = names(project.parameters), value = TRUE)
message("\tYAML config file contains ", length(recountIDs)," recount IDs.")


recountIDs.with.problems <- c("SRP003611" = "the number of samples drops after run merging")

## Optional: select a subset of the recountIDs
selectedRecountIDs <- c("SRP057196", "SRP042620")
# selectedRecountIDs <- setdiff(recountIDs, names(recountIDs.with.problems))

message("\tSelected ", length(selectedRecountIDs)," recount IDs: ", paste(collapse="; ", selectedRecountIDs))
#selectedRecountIDs <- recountIDs

# recountID <- "SRP042620"   ## Multi-group breast cancer
recountID <- "SRP057196"    # individual cells into all of the major neuronal, glial, and vascular cell types in the brain

# recountID <- "SRP003611"   # transcriptomes of 347 cells from 10 distinct populations in both of low-coverage (~0.27 million reads per cell) and high-coverage (~5 million reads per cell)
# recountID <- "SRP061240"   # several types from cancer (pancreatic, colorectal, prostat cancer) against Healthy control

#recountID <- "SRP006574"  # NOT working:  # this Project will constitute the classes that have IRanges object and then it will not have the same length for count Table such leads for STOP

# recountID <- "SRP062966"   # this project Boold disease
# recountID <- "SRP066834"   # cerebral cancer type.
# recountID <-  "SRP039694"  # hepatocellular carcinoma
# recountID <- "SRP008976"   # NOT working properly
######################### such these projects have promptly working #####################
# recountID <- "SRP042620"   ## Multi-group breast cancer. working excellent.
# recountID <- "SRP061240"   # several types from cancer (pancreatic, colorectal, prostat cancer) against Healthy control
# recountID <- "SRP062966"   # this project Blood disease. working excellent

# recountID <- "SRP056295"  # types of leukemia there is problem with random forest it say there is missing classes BUT actually there aren't any missing in the y class lable

# recountID <- "SRP003611"   # transcriptomes of 347 cells from 10 distinct populations in both of low-coverage (~0.27 million reads per cell) and high-coverage (~5 million reads per cell)
# recountID <- "SRP006574"   # NOT working:  # this Project will constitute the classes that have IRanges object and then it will not have the same length for count Table such leads for STOP
# recountID <- "SRP066834"   # Not working: cerebral cancer type.
# recountID <-  "SRP039694"  # Not working: hepatocellular carcinoma
# recountID <- "SRP008976"   # NOT working properly
# recountID <- "SRP006575"   # Not working


message("\tDefault recountID ", recountID)

################################################################
## Define general directories.
## Dataset-specific directories are defined later.


## END OF SCRIPT
#################################################################

