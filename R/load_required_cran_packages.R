
###################################################################
################### Required Library #############################
##################################################################

#' @title Load required packages from CRAN, after having installed them if necessary.
#' @author Mustafa AbuElQumsan and Jacques van Helden
#' @description Load a list of required RCRAN packages. For each specified package,
#' first check  if it is available on this R instance. If not, automatically install it.
#' @param packages a vector containing a list of required CRAN packages.
#' @return chack and then install all CRAN packages if it aren't installed.
#' @export
LoadRequiredCRANPackages <- function(packages, verbose = 0) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      if (verbose >= 1) { message("\tInstalling required CRAN library\t", pkg) }
      install.packages(pkg, dependencies = TRUE)
    }
    if (verbose >= 1) { message("\tLoading required CRAN library\t", pkg) }
    library(pkg, character.only = TRUE)
  }
}


