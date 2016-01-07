// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' Colsums of crossprod
//'
//' @description Fast computation of X %*% Y and colSums(X%*%Y)
//' @param X A matrix with dimensions k*n.
//' @param Y A matrix with dimenions n*m.
//' @param summit Logical. If \code{TRUE} return colsums of the resulting matrix.
//' @param transposeX Logical. If \code{TRUE} transpose X before multiplication.
//' @param transposeY Logical. If \code{TRUE} transpose Y before multiplication.
//' @return If \code{summit=0} a matrix with dimensions k * m else a vector of length m.
//' @author Thomas Alexander Gerds <tag@@biostat.ku.dk>
//' @export
// [[Rcpp::export]]
NumericMatrix cprod(NumericMatrix X, NumericMatrix Y, bool summit, bool transposeX, bool transposeY) {
  arma::mat A(X.begin(), X.nrow(), X.ncol(), false);
  arma::mat B(Y.begin(), Y.nrow(), Y.ncol(), false);
  arma::mat C;
  if (transposeX) 
    if (transposeY) 
      C = A.t() * B.t();
    else
      C = A.t() * B;
  else
    if (transposeY) 
      C = A * B.t();
    else
      C = A * B;
  int num_rows = C.n_rows;
  int num_cols = C.n_cols;
  if (summit){
    C = sum(C,0);
  }
  return wrap(C);
}
