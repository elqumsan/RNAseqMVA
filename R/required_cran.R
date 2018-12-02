
#require("roxygen2") ## MUSTAFA: do we really required roxygen2 and devtools to run the scripts ? They are required for building the package, but maybe not for running the analyses
#require("devtools") ## MUSTAFA: do we really required roxygen2 and devtools to run the scripts ? They are required for building the package, but maybe not for running the analyses

# loading required libraries
requiredCRAN <- c(
    'roxygen2',
    'dplyr', # Required for caret (depedency of dependency)
    'RCurl', # Undeclared dependency of recount
    'XML', # Undeclared dependency of recount
    'Rcpp',  # Wrongly declared dependency for xml2, itself required for roxygen2
    'roxygen2',
    'yaml',  # Required for the config file
    'devtools',
    'caret',
    'class',
    "broom",
    "roxygen2",
    "scatterplot3d",
    "e1071",
    "randomForest",
    "foreach",
    "doParallel",
    "pheatmap", # Heatmaps with extended control on dimensions
    "doMC")
# RequiredCRANPackages(requiredCRAN)

message("Required CRAN libraries:\n\t", paste(collapse="\n\t", requiredCRAN))

