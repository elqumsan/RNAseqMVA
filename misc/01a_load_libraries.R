### Loading the required external librares (CRAN, BioConductor) for RNAseqMVA ####

## Define the lists of required libraries
source('R/required_bioconductor.R')
source('R/required_cran.R')

## Read the functions to load the required libraries
source('R/load_required_cran_packages.R')
source('R/load_required_bioconductor_packages.R')

## Load required libraries
LoadRequiredCRANPackages(requiredCRAN)
LoadRequiredBioconductorPackages(requiredBioconductor)

## Load the RNAseqMVA package
require("RNAseqMVA")
