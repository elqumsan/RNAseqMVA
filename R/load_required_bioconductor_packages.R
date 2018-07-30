
###################################################################
################### Required Library #############################
##################################################################

#' @title Load required packages from BioConductor, afte having installed them if necessary.
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @description Load a list of required BioConductor packages. For each specified package,
#' first check  if it is available on this R instance. If not, automatically install it.
#' @param packages a vector containing a list of required BioConductor packages.
#' @return chack and then install all BioConductor packages if it aren't installed.
#' @export
LoadRequiredBioconductorPackages <- function(packages, verbose = 0) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      if (verbose >= 1) { message("\tInstalling required BioConductor library\t", pkg) }
      source('http://bioconductor.org/biocLite.R')
      # biocLite("BiocUpgrade") # note: I had to add this
      biocLite(pkg)
    }
    if (verbose >= 1) { message("\tLoading required BioConductor library\t", pkg) }
    library(pkg, character.only = TRUE)
  }
}

