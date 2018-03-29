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
selectedRecountIDs <- c("SRP042620", "SRP057196")
#selectedRecountIDs <- recountIDs

recountID <- "SRP042620"   ## Multi-group breast cancer
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


## Prefix for experiments with permuted class labels
## (negative controls to estimate random expectation)
perm.prefix <- parameters$perm_prefix


## Temporary, to get quick results
# parameters$classifiers <- c("knn")
# parameters$iterations <- 3


message.with.time <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d_%H%M%S"), "\t", ...)
}

################################################################
## Define directories

# Main directory should be adapted to the user's configuration
dir.main <- parameters$dir$main
tsv.dir <- paste(sep = "" ,parameters$dir$TSV,"/",recountID)

#classifier <- "knn"
## All other directories should be defined relative to dir.main
dir.scripts <- file.path(dir.main, "R")
dir.results <- file.path(parameters$dir$workspace, "results", parameters$recountID)
classifiers <- c("knn","rf", "svm")
dir.classifier <- file.path(dir.results, classifiers)

## Define the directories where tables and figures will be stored.
## one directory per classifer, with separate subdirectories for tables and figures.
classifier.dirs <- vector()
table.dirs <- vector()
figure.dirs <- vector()

detailFigures.dir <- vector()
detailTables.dir <- vector()

for (classifier in classifiers) {
  classifier.dirs[classifier] <- file.path(dir.results, classifier)
  dir.create(classifier.dirs[classifier], showWarnings = F, recursive = T)
  table.dirs[classifier] <- file.path(classifier.dirs[classifier], "tables")
  dir.create(table.dirs[classifier], showWarnings = F, recursive = T)

  detailTables.dir[classifier] <- file.path(table.dirs[classifier], "detailTables")
  dir.create(detailTables.dir[classifier],showWarnings = F, recursive = T)

  figure.dirs[classifier] <- file.path(classifier.dirs[classifier], "figures")
  dir.create(figure.dirs[classifier], showWarnings = F, recursive = T)

  detailFigures.dir[classifier] <- file.path(figure.dirs[classifier], "detailFigures")
  dir.create(detailFigures.dir[classifier] ,showWarnings = F, recursive = T)
}

## File to store a memory image
image.file <- file.path(dir.results, paste("RNA-seq_classifer_evaluation_", parameters$recountID, ".Rdata", sep = ""))

if (parameters$reload == TRUE) {
  ################################################################################
  ## Save an image of the memory, so I can reload it later to avoid re-running all the analyses
  parameters.current <- parameters # Keep current parameters to restore them after having loaded a memory image
  message.with.time("Loading memory image ")
  load(file = image.file)
  parameters <- parameters.current ## Reload current parameters (they might have been saved different in the memory image)
  rm(parameters.current)
}

## Load custom librarires
# source(file.path(dir.scripts, "misclassification_estimate.R") )
# source(file.path(dir.scripts, "load_counts.R"))
# source(file.path(dir.scripts, "required_libraries.R"))
# source(file.path(dir.scripts, "normalize_count_table.R"))
# source(file.path(dir.scripts, "deg_ordering.R"))
# source(file.path(dir.scripts, "one_experiment.R"))
# source(file.path(dir.scripts, "ErrorRateBoxPlot.R"))
# source(file.path(dir.scripts,"filterCountTable.R"))

# loading required libraries
requiredCRAN <- c('devtools', 'class', "randomForest","broom", "roxygen2","scatterplot3d","e1071")
RequiredCRANPackages(requiredCRAN)

## JvH: Mustafa, please add the other required packages, in particular recount
requiredBioconductor <- c("recount")
RequiredBioconductorPackages(requiredBioconductor)

################################################################
## TO CHECK LATER: DO wE STILL NEED THESE VARIABLES ???

## Directory for impact of Normalization and log2 into counts (and the study of its impact)
dir.NormImpact <- file.path(dir.results , paste("impact_of_normalisation_and_log2", sep = ""))
dir.create(dir.NormImpact, showWarnings = F, recursive = T)

## Directory for the visualization of Principal component for counts (and the study of its impact)
dir.visualisePCs <- file.path(dir.results , paste( "visualization_of_PCs", sep = ""))
dir.create(dir.visualisePCs, showWarnings = F, recursive = T)



## END OF SCRIPT
#################################################################

