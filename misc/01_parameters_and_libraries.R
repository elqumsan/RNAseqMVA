require("roxygen2")
require("devtools")
require("RNAseqMVA")
require("yaml")

## Read project-specific parameters from a yaml-formatted file.
## These parameters will then be used to overwrite the default parameters.
project.parameters <- yaml.load_file("misc/00_project_parameters.yml")

## Lad default parameters (must have been defined n the parameter field)
parameters <- project.parameters$default

## Get all recount IDs
recountIDs <- grep(pattern = "^SRP", x = names(project.parameters), value = TRUE)
message("YAML config file contains ", length(recountIDs)," recount IDs.")
selectedRecountIDs <- c("SRP042620", "SRP057196")
message("Selected ", length(selectedRecountIDs)," recount IDs: ", paste(collapse="; ", selectedRecountIDs))
#selectedRecountIDs <- recountIDs

recountID <- "SRP042620"   ## Multi-group breast cancer
message("Default recountID ", recountID)
# recountID <- "SRP057196"    # individual cells into all of the major neuronal, glial, and vascular cell types in the brain

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

# View(parameters)

#
#
# ## Parameters for knn
# parameters$knn = list(k=10) # the number of neighbours considered
#
# ## Parameters for svm
# parameters$svm = list(
#   kernel="linear",
#   scale=FALSE,
#   type="C-classification" ,
#   gamma = 1,
#   cost = 100)


## Temporary, to get quick results
# parameters$classifiers <- c("knn")
# parameters$iterations <- 3



# loading required libraries
requiredCRAN <- c('devtools', 'class', "randomForest","broom", "roxygen2","scatterplot3d","e1071")
RequiredCRANPackages(requiredCRAN)

## JvH: Mustafa, please add the other required packages, in particular recount
requiredBioconductor <- c("recount")
RequiredBioconductorPackages(requiredBioconductor)


## END OF SCRIPT
#################################################################

