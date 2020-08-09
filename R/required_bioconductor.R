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

message("Required Bioconductor libraries:\n\t", paste(collapse="\n\t", requiredBioconductor))

