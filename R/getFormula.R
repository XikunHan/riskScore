getFormula <- function(formula,object,verbose,...){
  if (missing(formula)){
    tryform <- try(formula <- eval(object[[1]]$call$formula),silent=TRUE)
    if ("try-error" %in% class(tryform))
        stop("Missing formula. Try to extract formula from (first element of) object failed.")
    if (match("formula",class(formula),nomatch=0)==0)
      stop("Argument formula is missing with no default.")
    else if (verbose>0)
      warning("Formula missing. Using formula from first model")
  }
  formula.names <- try(all.names(formula),silent=TRUE)
  if (!(formula.names[1]=="~")
      ||
      (match("$",formula.names,nomatch=0)+match("[",formula.names,nomatch=0)>0)){
    stop("Invalid specification of formula.\n Could be that you forgot the right hand side:\n ~covariate1 + covariate2 + ...?\nNote that any subsetting, ie data$var or data[,\"var\"], is not supported by this function.")
  }
  formula
}
