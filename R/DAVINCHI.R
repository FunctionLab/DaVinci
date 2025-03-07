#' DaVinci: Deconvolving lAtent Variables for Integrated Niche Cluster Identification
#' 
#' @useDynLib DaVinci, .registration=TRUE
NULL


#ArchR/R/GlobalDefaults.R
pkgs <- c("Rcpp",
"RcppArmadillo",
"RcppEigen",
"ggpubr",
"tripack",
"rsvd",
"cluster",
"mclust",
"dendextend",
"BiocNeighbors")


.onAttach <- function(libname, pkgname){
    #load packages
    packageStartupMessage("Loading Required Packages...")

    for (ii in seq_along(pkgs)){
         packageStartupMessage("\tLoading Package : ", pkgs[ii], " v", packageVersion(pkgs[ii]))
    tryCatch({
      suppressPackageStartupMessages(require(pkgs[ii], character.only=TRUE))
    }, error = function(e){
      packageStartupMessage("\tFailed To Load Package : ", pkgs[ii], " v", packageVersion(pkgs[ii]))
    })
    }#for ii

}#.onAttach


