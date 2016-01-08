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
NumericMatrix colSumsCrossprod(NumericMatrix X, NumericMatrix Y, bool transposeY){
  arma::mat A(X.begin(), X.nrow(), X.ncol(), false);
  arma::mat B(Y.begin(), Y.nrow(), Y.ncol(), false);
  arma::rowvec result;
  if (transposeY)
    result = arma::sum(A,1).t()*B.t();
  else
    result = arma::sum(A,1).t()*B;
  return wrap(result); 
}

//' @export
// [[Rcpp::export]]
NumericMatrix rowSumsCrossprod(NumericMatrix X, NumericMatrix Y, bool transposeY){
  arma::mat A(X.begin(), X.nrow(), X.ncol(), false);
  arma::mat B(Y.begin(), Y.nrow(), Y.ncol(), false);
  arma::rowvec result;
  if (transposeY)
      result = arma::sum(A,0)*B.t();
  else
      result = arma::sum(A,0)*B;
  return wrap(result); 
}


 
