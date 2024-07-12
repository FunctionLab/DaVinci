eigenMapMatMult <- function(A, B){
  .Call("cpp_eigenMapMatMult", A, B)
}#eigenMapMatMult


rcpp_sylvester <- function(left, right, total){
  .Call("cpp_rcpp_sylvester", left, right, total)
}#rcpp_sylvester


rcpp_shur <- function(x){
  .Call("cpp_rcpp_shur", x)
}#rcpp_shur


sylvester_pre <- function(left, Z2, T2, total){
  .Call("cpp_sylvester_pre", left, Z2, T2, total)
}#sylvester_pre


