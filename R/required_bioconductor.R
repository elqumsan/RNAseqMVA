
#require("roxygen2") ## MUSTAFA: do we really required roxygen2 and devtools to run the scripts ? They are required for building the package, but maybe not for running the analyses
#require("devtools") ## MUSTAFA: do we really required roxygen2 and devtools to run the scripts ? They are required for building the package, but maybe not for running the analyses

## Bioconductor packages required for the analyses.
## These cannot be declared with @import, since they are cannot be installed via install.packages()
requiredBioconductor <- c(
  "IRanges", # Seems to be required for derfinder but not automatically installed
  "recount",
  "vsn", ## For RNA-seq plots from Huber
  "DESeq2",
  "limma",
  "edgeR",
  "S4Vectors",
  "derfinder",
  "SummarizedExperiment")
# LoadRequiredBioconductorPackages(requiredBioconductor)

message("Required Bioconductor libraries:\n\t", paste(collapse="\n\t", requiredBioconductor))


## PROBLEM WITH derfinder
##  namespace ‘rlang’ 0.3.4 is already loaded, but >= 0.4.0 is required
## Error: package or namespace load failed for ‘derfinder’ in loadNamespace(i, c(lib.loc, .libPaths()), versionCheck = vI[[i]]):
## namespace ‘rlang’ 0.3.4 is already loaded, but >= 0.4.0 is required
