<p align="center">
  <img src="Logo.png" width="256">
</p>

# 

**DaVinci**: <ins>**D**</ins>econvolving l<ins>**a**</ins>tent <ins>**V**</ins>ariables for <ins>**i**</ins>ntegrated <ins>**n**</ins>iche <ins>**c**</ins>luster <ins>**i**</ins>dentification



## To install the R package from this github repository
 ```
 remotes::install_github("FunctionLab/DaVinci")
 ```
 ### If you have problems installing
 You can download the the repo to `working_directory` and use the following code block to inlcude all DaVinci functions in your working environment.
 ```
 library(Rcpp)
 library(mclust)
 
 script.path <- "working_directory/R/"
 script.list <- list.files(script.path)
 script.list <- setdiff(script.list, c("DAVINCHI.R", "import.R", "RcppExports.R"))
 script.list <- paste0(script.path, script.list)
 sapply(script.list, source)
 
 sourceCpp("working_directory/src/util.cpp")
 
 ```


## Tutorials

We provide step-by-step tutorials for applying DaVinci to a variety of spatial omics datasets. Detailed documentation is available [here](https://functionlab.github.io/DaVinci-docs/).
