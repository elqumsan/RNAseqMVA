message("Checking required libraries")

require("roxygen2") ## MUSTAFA: do we really required roxygen2 and devtools to run the scripts ? They are required for building the package, but maybe not for running the analyses
#require("devtools") ## MUSTAFA: do we really required roxygen2 and devtools to run the scripts ? They are required for building the package, but maybe not for running the analyses

## Bioconductor packages required for the analyses.
## These cannot be declared with @import, since they are cannot be installed via install.packages()
requiredBioconductor <- c(
  "IRanges", # Seems to be required for derfinder but not automatically installed
  "derfinder",
  "recount",
  "vsn", ## For RNA-seq plots from Huber
  "DESeq2",
  "limma",
  "edgeR",
  "S4Vectors",
  "SummarizedExperiment")
# RequiredBioconductorPackages(requiredBioconductor)
