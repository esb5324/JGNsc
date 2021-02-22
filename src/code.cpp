// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//' Multiply a number by two
//'
//' @param x A single integer.
//' @export
// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}

//' Thresholding a Partial correlation matrix
//'
//' @param x Gaussian matrix with samples x genes
//' @param y the estimated partial correlation matrix
//' @export
// [[Rcpp::export]]
double ThresholdPartCorr(NumericMatrix x, NumericMatrix y){
  int Ngene = y.ncol();
  double err =0;

  if (x.ncol() != y.ncol())
    Rf_error("Genes do not match.");
  for (int i=0; i<Ngene; i++){
    err += sum(pow((x(_,i) - y(_,i) * transpose(x) ),2));
  }
  return (err);
}

//' Thresholding a list of Partial correlation matrices
//'
//' @param x a list of Gaussian matrices with samples x genes
//' @param y a list of the estimated partial correlation matrices
//' @export
// [[Rcpp::export]]
double ThresholdPartCorrList(List x, List y){
  // x: the Gaussian matrix with samples x genes
  // y: the estimated partial correlation matrix
  arma::mat x0 = x[0];
  int Ngene = x0.n_cols;
  int K = x.length();
  double err =0;
  for(int k=0; k<K; k++){
    arma::mat xk = x[k];
    arma::mat yk = y[k];
    for (int i=0; i<Ngene; i++){
      err += sum(pow((xk.col(i) - xk * yk.col(i)),2));
    }
  }
  return (sqrt(err));
}

// //' rowSums cpp function
// //'
// //' @param X: matrix
// //' @export
// // [[Rcpp::export]]
// vec rowSums(const arma::mat & X){
//   int nRows = X.n_rows;
//   vec out(nRows);
//   for(int i = 0; i < nRows; i++){
//     out(i) = sum(X.row(i));
//   }
//   return(out);
// }
//
// //' colSums cpp function
// //'
// //' @param X: matrix
// //' @export
// // [[Rcpp::export]]
// vec colSums(const arma::mat & X){
//   int nCols = X.n_cols;
//   vec out(nCols);
//   for(int i = 0; i < nCols; i++){
//     out(i) = sum(X.col(i));
//   }
//   return(out);
// }
