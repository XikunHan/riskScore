#' Probability Predictions
#' 
#' Function to extract probabilistic event status predictions from
#' various diagnostic and prognostic models with binary status response. The
#' function has a speficic method depending on the 'class' of the object.
#' 
#' The function delivers predicted probabilities tailored for the model
#' performance measures of the package. These probabilities are extracted from
#' a fitted model of class \code{CLASS} with the function
#' \code{predictStatusProb.CLASS}. See \code{help(Roc)} for details.
#' 
#' @aliases predictStatusProb predictStatusProb.randomForest
#' predictStatusProb.lrm predictStatusProb.default predictStatusProb.glm
#' predictStatusProb.rpart predictStatusProb.rfsrc predictStatusProb.numeric 
#' @usage
#' \method{predictStatusProb}{glm}(object,newdata,...)
#' @param object A model from which predicted probabilities can be
#' extracted for the indiviuals in newdata.
#' @param newdata A data frame containing data for which the \code{object}
#' can provide predict probabilities. In medical
#' applications \code{newdata} will typically consist of the data of 
#' patients whose data were not used for building the model.
#' @param ... Additional arguments that are passed on to the current
#' method.
#' @return A vector with the predicted status probability for each row in
#' \code{NROW(newdata)}.
#' 
#' @note It is rather easy to write a new predictStatusProb method, see \code{help(Roc)}. However,
#' if you do not succeed, please send me an email.
#' 
#' The performance, in particular when doing cross-validation where the model
#' is evaluated many times, can be improved by supressing in the call to the
#' model all the computations that are not needed for probability prediction,
#' for example standard error calculations.
#'
#' 
#' @usage predictStatusProb(object,newdata, ...)
#' @examples
#' library(rms)
#' set.seed(7)
#' x <- abs(rnorm(20))
#' d <- data.frame(y=rbinom(20,1,x/max(x)),x=x,z=rnorm(20))
#' nd <- data.frame(y=rbinom(8,1,x/max(x)),x=abs(rnorm(8)),z=rnorm(8))
#' fit <- lrm(y~x+z,d)
#' predictStatusProb(fit,newdata=nd)
#' 
#' @author Thomas A. Gerds \email{tag@@biostat.ku.dk}
#' @seealso \code{\link{predict}},\code{\link{Roc}}
#' @keywords models
#' @export predictStatusProb
predictStatusProb <- function(object,...){
  UseMethod("predictStatusProb")
}

##' @export
predictStatusProb.numeric <- function(object,newdata,...){
  stopifnot(NROW(object)==NROW(newdata))
  p <- object
  class(p) <- "predictStatusProb"
  p
}

##' @export
predictStatusProb.formula <- function(object,newdata,...){
    ff <- update.formula(object,"NULL~.")
    if (length(all.vars(ff))==1){
        p <- stats::model.frame(ff,newdata)[[1]]
        class(p) <- "predictStatusProb"
        p
    } else{
        fit <- glm(object,data=newdata,family="binomial")
        predictStatusProb(fit,newdata=newdata,...)
    }
}

##' @export
predictStatusProb.double <- function(object,newdata,...){
    stopifnot(NROW(object)==NROW(newdata))
    p <- object
    class(p) <- "predictStatusProb"
    p
}

## predictStatusProb.NULL <- function(object,newdata,...){
## runif(NROW(newdata))
## }
##' @export
predictStatusProb.glm <- function(object,newdata,...){
    if (object$family$family=="binomial")
        p <- as.numeric(stats::predict(object,newdata=newdata,type="response"))
    else{ stop("Currently only the binomial family is implemented for predicting a status from a glm object.")
      }
    class(p) <- "predictStatusProb"
    p
}
##' @export
predictStatusProb.BinaryTree <- function(object,newdata,...){
    treeresponse <- party::treeresponse
    p <- sapply(treeresponse(object,newdata=newdata),function(x)x[1])
    class(p) <- "predictStatusProb"
    p
}

##' @export
predictStatusProb.lrm <- function(object,newdata,...){
  p <- as.numeric(stats::predict(object,newdata=newdata,type="fitted"))
  class(p) <- "predictStatusProb"
  p
}

##' @export
predictStatusProb.rpart <- function(object,newdata,...){
  p <- as.numeric(stats::predict(object,newdata=newdata,type="prob")[,2,drop=TRUE])
  class(p) <- "predictStatusProb"
  p
}

##' @export
predictStatusProb.randomForest <- function(object,newdata,...){
  stopifnot(!missing(newdata))
  p <- as.numeric(stats::predict(object,newdata=newdata,type="prob")[,2,drop=TRUE])
  class(p) <- "predictStatusProb"
  p
}

##' @export
predictStatusProb.rfsrc <- function(object, newdata, ...){
  p <- as.numeric(stats::predict(object,newdata=newdata,importance="none",...)$predicted[,2])
  class(p) <- "predictStatusProb"
  p
}



