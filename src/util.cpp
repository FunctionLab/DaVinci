// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace arma;



//' Fast Matrix multiplication
//' @param A 1st matrix
//' @param B 2nd matrix
//' @export
// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B;

  return Rcpp::wrap(C);
}



arma::mat rcpp_sylvester(arma::mat left, arma::mat right, arma::mat total){
  return syl(left, right, -total);
}//rcpp_sylvester


//' Schur decomposition
//' @param x input matrix
//' @export
// [[Rcpp::export]]
Rcpp::List rcpp_shur(arma::mat x){
  arma::mat U;
  arma::mat S;
  bool info = schur(U, S, x);

  return Rcpp::List::create(Rcpp::Named("U") = U, Rcpp::Named("S") = S, Rcpp::Named("info") = info);
}//rcpp_shur



//' Solve the sylvester equation
//' @param left
//' @param Z2
//' @param T2
//' @param total
//' @export
// [[Rcpp::export]]
arma::mat sylvester_pre(arma::mat left, arma::mat Z2, arma::mat T2, arma::mat total){
  arma::mat Z1, T1;
  bool info1 = schur(Z1, T1, left);
  //bool info2 = schur(Z2, T2, right);

  char trana = 'N';
  char tranb = 'N';
  int isgn = +1;
  
  blas_int m = (blas_int)T1.n_rows;
  blas_int n = (blas_int)T2.n_cols;
  blas_int lda = m;
  blas_int ldb = n;
  blas_int ldc = m;

  //Rcout << "0" << "\n";

  double scale = 0.0;
  int info = 0;

  arma::mat Y = Z1.t()*total*Z2;

  //LAPACKE_dtrsyl(info, trana, tranb, isgn, m, n, T1.memptr(), m, T2.memptr(), n, Y.memptr(), m, scale);
  //F77_CALL(dtrsyl)(&trana, &tranb, &isgn, &m, &n, T1.memptr(), &m, T2.memptr(), &n, Y.memptr(), &m, &scale, &info);
  //dtrsyl_(&trana, &tranb, &isgn, &m, &n, T1.memptr(), &m, T2.memptr(), &n, Y.memptr(), &m, &scale, &info);

  arma::lapack::trsyl(&trana, &tranb, &isgn, &m, &n, 
                       T1.memptr(), &lda, 
                       T2.memptr(), &ldb, 
                       Y.memptr(), &ldc, 
                       &scale, &info);
                       
  //Y /= (-scale);
  Y /= (scale);

  arma::mat X = Z1*Y*Z2.t();
  return X;
}//sylvester_pre




