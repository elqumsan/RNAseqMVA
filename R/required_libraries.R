
###################################################################
################### Required Library #############################
##################################################################

#' @title Load required packages from CRAN, and install them if necessary. 
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @description Load a list of required RCRAN packages. For each specified package, 
#' first check  if it is available on this R instance. If not, automatically install it. 
#' @param packages a vector containing a list of required CRAN packages.
#' @export
RequiredCRANPackages <- function (packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

#' @title Load required packages from BioConductor, and install them if necessary. 
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @description Load a list of required BioConductor packages. For each specified package, 
#' first check  if it is available on this R instance. If not, automatically install it. 
#' @param packages a vector containing a list of required BioConductor packages.
#' @export
RequiredBioconductorPackages <-function (packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      source('http://bioconductor.org/biocLite.R')
      # biocLite("BiocUpgrade") # note: I had to add this
      biocLite(pkg)
    }
    library(pkg, character.only = TRUE)
  }
} 

