// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include </mnt/sw/nix/store/mm1sfzml36dranq9gki8h7lmcfiy4bgz-r-4.1.3-view/rlib/R/library/RcppEigen/include/Eigen/src/misc/lapacke.h>
using namespace Rcpp;
using namespace arma;



// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B;
  
  return Rcpp::wrap(C);
}


// [[Rcpp::export]]
arma::mat rcpp_sylvester(arma::mat left, arma::mat right, arma::mat total){
  return syl(left, right, -total);
}//rcpp_sylvester


// [[Rcpp::export]]
Rcpp::List rcpp_shur(arma::mat x){
  arma::mat U;
  arma::mat S;
  bool info = schur(U, S, x);
  
  return Rcpp::List::create(Rcpp::Named("U") = U, Rcpp::Named("S") = S, Rcpp::Named("info") = info);
}//rcpp_shur



// [[Rcpp::export]]
arma::mat sylvester_pre(arma::mat left, arma::mat Z2, arma::mat T2, arma::mat total){
  arma::mat Z1, T1;
  bool info1 = schur(Z1, T1, left);
  //bool info2 = schur(Z2, T2, right);
  
  char trana = 'N';
  char tranb = 'N';
  int isgn = +1;
  int m = T1.n_rows;
  int n = T2.n_cols;
  
  //Rcout << "0" << "\n";
  
  double scale = 0.0;
  int info = 0;
  
  arma::mat Y = Z1.t()*total*Z2;
  
  //LAPACKE_dtrsyl(info, trana, tranb, isgn, m, n, T1.memptr(), m, T2.memptr(), n, Y.memptr(), m, scale);
  F77_CALL(dtrsyl)(&trana, &tranb, &isgn, &m, &n, T1.memptr(), &m, T2.memptr(), &n, Y.memptr(), &m, &scale, &info);
  
  //Y /= (-scale);
  Y /= (scale);
  
  arma::mat X = Z1*Y*Z2.t();
  return X;
}//sylvester_pre




