// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// colSumsCrossprod
NumericMatrix colSumsCrossprod(NumericMatrix X, NumericMatrix Y, bool transposeY);
RcppExport SEXP riskScore_colSumsCrossprod(SEXP XSEXP, SEXP YSEXP, SEXP transposeYSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< bool >::type transposeY(transposeYSEXP);
    __result = Rcpp::wrap(colSumsCrossprod(X, Y, transposeY));
    return __result;
END_RCPP
}
// rowSumsCrossprod
NumericMatrix rowSumsCrossprod(NumericMatrix X, NumericMatrix Y, bool transposeY);
RcppExport SEXP riskScore_rowSumsCrossprod(SEXP XSEXP, SEXP YSEXP, SEXP transposeYSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< bool >::type transposeY(transposeYSEXP);
    __result = Rcpp::wrap(rowSumsCrossprod(X, Y, transposeY));
    return __result;
END_RCPP
}
