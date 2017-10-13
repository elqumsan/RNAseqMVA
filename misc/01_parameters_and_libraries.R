## Load the parameters and libraries that will be used for the different other scripts.
recountIDs <- c("SRP042620", ## Multi-group breast cancer
                "SRP057196",
                "SRP056295", ## 525 runs, human leukemia
                # crash                "SRP042161", ## 875 runs, 430 single glioblastoma cells isolated from 5 individual tumors and 102 single cells from gliomasphere cells lines generated using SMART-seq
                "SRP041736", # transcriptomes of 347 cells from 10 distinct populations in both of low-coverage (~0.27 million reads per cell) and high-coverage (~5 million reads per cell)
                "SRP035988",  ## Psoriasis
                "ERP003613" # Tissue samples from 95 human individuals
)




<<<<<<< HEAD

=======
>>>>>>> 78b4dbb65c61735789ecc021a414a85aadfbbfa6
recountID <- "SRP042620" ## Multi-group breast cancer
#recountID <- "SRP041736" # transcriptomes of 347 cells from 10 distinct populations in both of low-coverage (~0.27 million reads per cell) and high-coverage (~5 million reads per cell)


## Running parameters
parameters <- list(
  recountID = recountID,
  classColumn = "tissue",
  mergeRuns = TRUE, ## Whether or not to merge runs per sample
  sampleIdColumn = "geo_accession",
  minSamplesPerClass = 10,
  #  iterations = 50, ## Number of iterations for the classiifers
  iterations = 10, ## Number of iterations for the classiifers
  reload = FALSE,
  compute = TRUE, ## If FALSE, do not run the heavy computations, just generate the pictures and save tables
  save.tables = TRUE, ## If TRUE, save tab-delimited files with the results
  save.image = TRUE,
  #  classifiers = c("lda"),
  classifiers = c("knn", "rf", "lda","svm"),
  data.types = c("raw", "norm", "log2"),
  deg.methods = c("DESeq2", "edgeR", "randomized"),
  ## Note: for the number of variables I choose a regular spacing around the total number of samples in order to detect overfitting effects
  nb.variables = c( 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 300, 400, 500, 1000, 2000, 5000, 10000), ## Number of variables for variable ordering test
  trainingProportion = 2/3, # ratio of spliting the data set to training and testing sets
  permute = c(FALSE, TRUE),
  verbose = TRUE
)


## Parameters for knn
parameters$knn = list(k=10) # the number of neighbours considered

## Parameters for svm
parameters$svm = list(
  kernel="linear",
  scale=FALSE,
  type="C-classification" ,
  gamma = 1,
  cost = 100)


if (parameters$recountID == "SRP042620") {
  ## Psoriasis dataset
  parameters$classColumn <- "tissue"
} else if (parameters$recountID == "SRP041736") {
  parameters$sampleIdColumn <- "sample"
  parameters$classColumn <- NA
  stop("We cannot analyse this dataset because the pheno table does not contain any info about the sample classes")

} else if (parameters$recountID == "SRP057196") {
  ## Running parameters for SRP057196
  ## NOTE FOR MUSTAFA: For this dataset, the groups are defined
  ## based on a combination of two columns of the pheno table:
  ## "tissue" and "cell.type".
  ## I improved filterCountTable() in order to support combinations of columns.
  ## If several columns are specified in classColumn, the classes
  ## are built by pasting the specified columns.
  parameters$classColumn <- c("tissue", "cell.type")
} else if (parameters$recountID == "SRP035988") {
  parameters$classColumn <- "tissue.type"
}


## Prefix for experiments with permuted class labels
## (negative controls to estimate random expectation)
perm.prefix <- "permLabels"


## Temporary, to get quick results
# parameters$classifiers <- c("knn")
# parameters$iterations <- 3


message.with.time <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d_%H%M%S"), "\t", ...)
}

################################################################
## Define directories

## Main directory should be adapted to the user's configuration
dir.main <- "~/RNAseqMVA"

#classifier <- "knn"
## All other directories should be defined relative to dir.main
dir.scripts <- file.path(dir.main, "R")
dir.results <- file.path(dir.main, "results", parameters$recountID)
classifiers <- c("knn","rf", "svm")
dir.classifier <- file.path(dir.results, classifiers)

## Define the directories where tables and figures will be stored.
## one directory per classifer, with separate subdirectories for tables and figures.
classifier.dirs <- vector()
table.dirs <- vector()
figure.dirs <- vector()
for (classifier in classifiers) {
  classifier.dirs[classifier] <- file.path(dir.results, classifier)
  dir.create(classifier.dirs[classifier], showWarnings = F, recursive = T)
  table.dirs[classifier] <- file.path(classifier.dirs[classifier], "tables")
  dir.create(table.dirs[classifier], showWarnings = F, recursive = T)
  figure.dirs[classifier] <- file.path(classifier.dirs[classifier], "figures")
  dir.create(figure.dirs[classifier], showWarnings = F, recursive = T)
}

## File to store a memory image
image.file <- file.path(dir.results, paste("RNA-seq_classifer_evaluation_", parameters$recountID, ".RData", sep = ""))

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
source(file.path(dir.scripts, "misclassification_estimate.R") )
source(file.path(dir.scripts, "load_counts.R"))
source(file.path(dir.scripts, "required_libraries.R"))
source(file.path(dir.scripts, "normalize_count_table.R"))
source(file.path(dir.scripts, "deg_ordering.R"))
source(file.path(dir.scripts, "one_experiment.R"))
source(file.path(dir.scripts, "ErrorRateBoxPlot.R"))
source(file.path(dir.scripts,"filterCountTable.R"))

# loading required libraries
requiredCRAN <- c('class', "randomForest","broom", "devtools")
RequiredCRANPackages(requiredCRAN)

## JvH: Mustafa, please add the other required packages, in particular recount
requiredBioconductor <- c("recount")
RequiredBioconductorPackages(requiredBioconductor)

################################################################
## TO CHECK LATER: DO wE STILL NEED THESE VARIABLES ???

## Directory for the normalized counts (and the study of its impact)
dir.NormImpact <- file.path(dir.results , paste(parameters$recountID, "_impacte_of_normalisation", sep = ""))
dir.create(dir.NormImpact, showWarnings = F, recursive = T)

## Directory for the log2-transformed counts (and the study of its impact)
dir.log2Impact <- file.path(dir.results , paste(parameters$recountID, "_impacte_of_log2", sep = ""))
dir.create(dir.log2Impact, showWarnings = F, recursive = T)

## END OF SCRIPT
#################################################################

