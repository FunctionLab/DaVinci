// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <R_ext/Lapack.h>
using namespace Rcpp;
using namespace arma;




extern "C" void F77_NAME(dtrsyl)(const char*, const char*, const int*, const int*, const int*,
                        const double*, const int*, const double*, const int*,
                        double*, const int*, double*, int*);


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
  int m = T1.n_rows;
  int n = T2.n_cols;

  //Rcout << "0" << "\n";

  double scale = 0.0;
  int info = 0;

  arma::mat Y = Z1.t()*total*Z2;

  //LAPACKE_dtrsyl(info, trana, tranb, isgn, m, n, T1.memptr(), m, T2.memptr(), n, Y.memptr(), m, scale);
  //F77_CALL(dtrsyl)(&trana, &tranb, &isgn, &m, &n, T1.memptr(), &m, T2.memptr(), &n, Y.memptr(), &m, &scale, &info);

  dtrsyl_(&trana, &tranb, &isgn, &m, &n, T1.memptr(), &m, T2.memptr(), &n, Y.memptr(), &m, &scale, &info);

  //Y /= (-scale);
  Y /= (scale);

  arma::mat X = Z1*Y*Z2.t();
  return X;
}//sylvester_pre




