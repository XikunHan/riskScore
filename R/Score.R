#' @importFrom foreach "%dopar%" foreach
#' @importFrom survival Surv
#' @importFrom prodlim Hist jackknife prodlim sindex
#' @importFrom data.table data.table :=
#' @importFrom grDevices col2rgb gray
#' @importFrom graphics abline axis box legend lines mtext par plot points segments text title
#' @importFrom stats model.frame model.response as.formula coef family formula median model.matrix na.fail na.omit pnorm predict quantile rbinom rexp runif sd smooth terms time update update.formula var wilcox.test qnorm
#' @importFrom utils capture.output head select.list
##' @export
Score <- function(object,...){
  UseMethod("Score",object=object)
}
##' List method to Score risk predictions
##'
##' Describe details
##' @title Score risk predictions
##' @aliases Score
##' @param object List of risk predictions (see details and examples). 
##' @param formula A formula which identifies the outcome (left hand side). For right censored outcome, 
##' the right hand side of the formula is used to estimate the IPCW model. 
##' @param data Data set or table in which the formula can be interpreted.
##' @param metrics Character vector specifying which metrics to apply. Implemented are \code{"auc"} and \code{"Brier"}.
##' @param plots Character vector specifying which plots to prepare. 
##' @param cause Event of interest. Used for binary outcome \code{Y} to specify that risks are risks of the event \code{Y=event}
##' and for competing risks outcome to specify the cause of interest.
##' @param times For survival and competing risks outcome: list of prediction horizons. All times which are greater
##' than the maximal observed time in the data set are removed.
##' @param landmarks Not yet implemented.
##' @param useEventTimes If \code{TRUE} add all unique event times to argument \code{times}.
##' @param nullModel If \code{TRUE} add the null model which ignores the covariates and predicts the prevalence for all subjects. 
##' @param test If \code{TRUE} compute confidence intervals. Also do model comparisons specified in \code{dolist}.
##' @param alpha Level of significance.
##' @param dolist Vector of integers specifying which risks to compare to all other risks in \code{object}.
##' @param censMethod Method for dealing with right censored data. Either \code{"ipcw"} or \code{"pseudo"}.
##' @param censModel Model for estimating inverse probability of censored weights.
##' @param splitMethod Method for cross-validation.
##' @param B Number of cross-validation steps.
##' @param M Size of subsamples for cross-validation. If specified it has to be an integer smaller than the size of \code{data}.
##' @param trainseeds Seeds for training models during cross-validation.
##' @param ... Not yet used.
##' @return Data table with scores and tests.
##' @examples
##' # binary outcome
##' library(lava)
##' set.seed(18)
##' learndat <- sampleData(100,outcome="binary")
##' testdat <- sampleData(40,outcome="binary")
##'
##' # score logistic regression models
##' lr1 = glm(Y~X1+X2+X7+X9,data=learndat,family=binomial)
##' lr2 = glm(Y~X3+X5+X6,data=learndat,family=binomial)
##' Score(list("LR(X1+X2+X7+X9)"=lr1,"LR(X3+X5+X6)"=lr2),formula=Y~1,data=testdat)
##' 
##' # compute AUC for a list of continuous markers
##' markers = as.list(testdat[,1:5])
##' u=Score(markers,formula=Y~1,data=testdat,metrics=c("auc"))
##'
##' # cross-validation
##' lr1a = glm(Y~X6,data=learndat,family=binomial)
##' lr2a = glm(Y~X7+X8+X9,data=learndat,family=binomial)
##' Score(list("LR1"=lr1a,"LR2"=lr2a),formula=Y~1,data=learndat,splitMethod="bootcv",B=3)
##'
##' # survival outcome
##' 
##' # Score Cox regression models
##' library(survival)
##' library(rms)
##' library(prodlim)
##' set.seed(18)
##' trainSurv <- sampleData(100,outcome="survival")
##' testSurv <- sampleData(40,outcome="survival")
##' cox1 = coxph(Surv(time,event)~X1+X2+X7+X9,data=trainSurv)
##' cox2 = coxph(Surv(time,event)~X3+X5+X6,data=trainSurv)
##' Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
##'       formula=Surv(time,event)~1,data=testSurv,test=FALSE,times=c(5,8))
##'
##' # time-dependent AUC for list of markers
##' survmarkers = as.list(testSurv[,1:5])
##' Score(survmarkers,
##'       formula=Surv(time,event)~1,metrics="auc",data=testSurv,
##' test=TRUE,times=c(5,8))
##' 
##' # compare models on test data
##' Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
##'       formula=Surv(time,event)~1,data=testSurv,test=TRUE,times=c(5,8))
##'
##' # crossvalidation models in traindata
##' Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
##'       formula=Surv(time,event)~1,data=trainSurv,test=TRUE,times=c(5,8),
##' splitMethod="bootcv",B=3)
##'
##' # restrict number of comparisons
##' Score(list("Cox(X1+X2+X7+X9)"=cox1,"Cox(X3+X5+X6)"=cox2),
##'       formula=Surv(time,event)~1,data=trainSurv,dolist=2,
##' nullModel=FALSE,test=TRUE,times=c(5,8),splitMethod="bootcv",B=3)
##'
##' # competing risks outcome
##' set.seed(18)
##' trainCR <- sampleData(40,outcome="competing.risks")
##' testCR <- sampleData(40,outcome="competing.risks")
##' library(riskRegression)
##' library(cmprsk)
##' # Cause-specific Cox regression
##' csc1 = CSC(Hist(time,event)~X1+X2+X7+X9,data=trainCR)
##' csc2 = CSC(Hist(time,event)~X3+X5+X6,data=trainCR)
##' # Fine-Gray regression
##' fgr1 = FGR(Hist(time,event)~X1+X2+X7+X9,data=trainCR,cause=1)
##' fgr2 = FGR(Hist(time,event)~X3+X5+X6,data=trainCR,cause=1)
##' Score(list("CSC(X1+X2+X7+X9)"=csc1,"CSC(X3+X5+X6)"=csc2,
##'            "FGR(X1+X2+X7+X9)"=fgr1,"FGR(X3+X5+X6)"=fgr2),
##'       formula=Hist(time,event)~1,data=testCR,test=TRUE,times=c(5,8))
##' 
##' @author Thomas A Gerds \email{tag@@biostat.ku.dk} and Paul Blanche \email{paul.blanche@@univ-ubs.fr}
##' @export 
Score.list <- function(object,
                       formula,
                       data,
                       metrics=c("auc","brier"),
                       plots=c("roc","calibration"),
                       cause=1,
                       times,
                       landmarks,
                       useEventTimes=FALSE,
                       nullModel=TRUE,
                       test=TRUE,
                       alpha=0.05,
                       dolist,
                       censMethod="ipcw",
                       censModel="cox",
                       splitMethod,
                       B,
                       M,
                       trainseeds,
                       ...){

    id=time=status=id=WTi=b=time=status=model1=model2=p=model=NULL

    # -----------------parse arguments and prepare data---------
    # {{{ Response
    if (missing(formula)){stop("Argument formula is missing.")}    
    formula.names <- try(all.names(formula),silent=TRUE)
    if (!(formula.names[1]=="~")
        ||
        (match("$",formula.names,nomatch=0)+match("[",formula.names,nomatch=0)>0)){
        stop("Invalid specification of formula.\n Could be that you forgot the right hand side:\n ~covariate1 + covariate2 + ...?\nNote that any subsetting, ie data$var or data[,\"var\"], is not supported by this function.")
    }
    if (missing(data)){stop("Argument data is missing.")}
    data <- data.table(data)
    responseFormula <- stats::update(formula,~1)
    ## if (missing(event)) event <- 1
    responsevars <- all.vars(responseFormula)
    response <- getResponse(formula=responseFormula,
                            data=data,
                            cause=cause,
                            vars=responsevars)
    response.dim <- NCOL(response)
    responseType <- attr(response,"model")
    # add null model and find names for the object
    if (nullModel==TRUE){
        nullobject <- getNullModel(formula=formula,data=data,responseType=responseType)
    } else{
          nullobject <- NULL
      }
    ## put ReSpOnSe for binary and (time, event, status) in the first column(s) 
    ## data[,eval(responsevars):=NULL]
    data <- cbind(response,data)
    if (responseType=="survival")
        formula <- stats::update(formula,"Hist(time,status)~.")
    if (responseType=="competing.risks")
        formula <- stats::update(formula,"Hist(time,event)~.")
    N <- NROW(response)
    predictHandlerFun <- switch(responseType,
                                "survival"="predictSurvProb",
                                "competing.risks"="predictEventProb",
                                "binary"="predictStatusProb",
                                stop("Dont know how to predict response of type ",responseType))
    censType <- attr(response,"cens.type")
    if (is.null(censType)) censType <- "uncensoredData"
    # }}}
    # {{{ SplitMethod
    splitMethod <- getSplitMethod(splitMethod=splitMethod,B=B,N=N,M=M)
    B <- splitMethod$B
    splitIndex <- splitMethod$index
    do.resample <- !(is.null(splitIndex))
    # }}}
    # {{{ Checking the models
    # for predictHandlerFunction
    allmethods <- utils::methods(predictHandlerFun)
    ## wantedMethods <- lapply(object,function(o){
    ## candidateMethods <- paste(predictHandlerFun,class(o),sep=".")
    ## if (all(match(candidateMethods,allmethods,nomatch=0)==0))
    ## stop(paste("Could not find ",predictHandlerFun," method for ",paste(class(o),collapse=" ,"),sep=""))
    ## })
    # checking the models for compatibility with resampling
    if (is.null(names(object))){
        names(object) <- sapply(object,function(o)class(o)[1])}
    else {
        names(object)[(names(object)=="")] <- sapply(object[(names(object)=="")],function(o)class(o)[1])
    }
    names.object <- names(object) <- make.unique(names(object))
    NF <- length(object)
    if (!is.null(nullobject)) {
        mlevs <- 0:NF
        mlabels <- c(names(nullobject),names(object))
    } else{
          mlevs <- 1:NF
          mlabels <- names(object)
      }
    if (do.resample){
        nix <- lapply(1:length(object),function(f){
                                  fit <- object[[f]]
                                  if(is.null(fit$call))
                                      stop(paste("model",names(object)[f],"does not have a call argument."))
                                  ## else fit$call$data <- NULL
                                  ## fit
                              })
        ## names(object) <- names.object
    }
    if ((NF+length(nullobject))<=1) dolist <- NULL 
    if (test==FALSE) dolist <- NULL
    ## test <- FALSE
    if (test==TRUE && missing(dolist)){
        if (is.null(nullobject)) {
            dolist <- NF:2
        } else{
              dolist <- NF:1
          }
    }
    # }}}
    # {{{ add id *before* ordering the data. data have to be ordered for pec::ipcw
    data[,id:=1:N]
    # }}}
    # {{{ Evaluation landmarks and horizons (times)
    if (responseType %in% c("survival","competing.risks")){
        ## in case of a tie, events are earlier than right censored
        data.table::setorder(data,time,-status)
        eventTimes <- unique(data[,time])
        maxtime <- eventTimes[length(eventTimes)]
        if (missing(landmarks)){
            start <- 0
            if (missing(times)){
                if (useEventTimes==TRUE)
                    times <- unique(c(start,eventTimes))
                else{
                    ## times <- seq(start,maxtime,(maxtime - start)/100)
                    times <- median(eventTimes)
                }
            } else{
                  if (useEventTimes==TRUE) 
                      times <- sort(c(start,unique(times),eventTimes))
                  else
                      ## times <- sort(unique(c(start,times)))
                      times <- sort(unique(times))
              }
            stopifnot(sum(times<=maxtime)>0)
            times <- times[times<=maxtime]
            NT <-  length(times)
        }
        else{
            stop("Landmark updating not yet implemented.")
        }
    }
    else{
        if (!missing(times)) warning("Function 'Score': Response type is not time-to-event: argument 'times' will be ignored.",call.=FALSE)
        times <- NULL
        NT <- 1
    }
    # }}}
    # ----------------------------find metrics and plots ----------------------
    # {{{
    
    ## Metrics <- lapply(metrics,grep,c("AUC","Brier"),ignore.case=TRUE,value=TRUE)
    metrics[grep("^auc$",metrics,ignore.case=TRUE)] <- "AUC"
    metrics[grep("^brier$",metrics,ignore.case=TRUE)] <- "Brier"
    Plots <- lapply(plots,grep,c("Roc","Cal"),ignore.case=TRUE,value=TRUE)
    # }}}
    # -----------------IPCW outside loop -----------------------
    # {{{ 
    if (responseType %in% c("survival","competing.risks")){
        if (censType=="rightCensored"){
            if ("outside" %in% censMethod){
                if ("ipcw" %in% censMethod){
                    Weights <- getCensoringWeights(formula=formula,
                                                   data=testdata,
                                                   response=response,
                                                   times=times,
                                                   censModel=censModel,
                                                   responseType=responseType)
                }
                else{
                    censMethod <- "jackknife.pseudo.values"
                    pseudoResponse <- getPseudoValues(formula=formula,data=testdata,responseType=responseType,times=times,cause=cause)
                }
            }else{
                 censMethod <- c(censMethod,"inside")
             }
        }
        else{
            if (censType=="uncensored"){
                censMethod <- c("none","inside")
                Weights <- NULL
            }
            else{
                stop("Cannot handle this type of censoring.")
            }
        }
    }
    # }}}
    # -----------------performance program----------------------
    # {{{ 
    trainModel <- function(model,data){
        model$call$data <- data
        try(eval(model$call),silent=TRUE)
    }
    computePerformance <- function(object,
                                   nullobject,
                                   testdata,
                                   traindata=NULL,
                                   trainseed=NULL,
                                   metrics,
                                   times,
                                   cause,
                                   test,
                                   alpha,
                                   dolist,
                                   NF,
                                   NT){
        N <- NROW(testdata)
        # split into response and predictors
        response <- testdata[,1:response.dim,with=FALSE]
        response[,id:=testdata[["id"]]]
        X <- testdata[,-c(1:response.dim),with=FALSE]
        # {{{ -----------------IPCW inner loop-----------------------
        if (responseType %in% c("survival","competing.risks")){
            if ("inside" %in% censMethod){
                if (censType=="rightCensored"){
                    if ("ipcw" %in% censMethod){
                        Weights <- getCensoringWeights(formula=formula,
                                                       data=testdata,
                                                       response=response,
                                                       times=times,
                                                       censModel=censModel,
                                                       responseType=responseType)
                        ## add subject specific weights here, and time specific weights later
                        response[,WTi:=Weights$IPCW.subjectTimes]
                    }else{
                         censMethod <- "jackknife.pseudo.values"
                         pseudoResponse <- getPseudoValues(formula=formula,data=testdata,responseType=responseType,times=times,cause=cause)
                     }
                } else{
                      if (censType=="uncensored"){
                          Weights <- list(IPCW.times=rep(1,NT),
                                          IPCW.subjectTimes=matrix(1,ncol=NT,nrow=N))
                          Weights$method <- "marginal"
                          response[,WTi:=1]
                      } else{
                            stop("Cannot handle this type of censoring.")
                        }
                  }
            }
        } else{
              Weights <- NULL
          }
        # }}}
        # extract predictions as data.table
        args <- switch(predictHandlerFun,"predictStatusProb"={list(newdata=X)},
                       "predictSurvProb"={list(newdata=X,times=times)},
                       "predictEventProb"={list(newdata=X,times=times,cause=cause)},
                       stop("Unknown predictHandler."))
        pred <- data.table::rbindlist(lapply(mlevs, function(f){
                                                 if (f!=0 && any(c("integer","factor","numeric","matrix") %in% class(object[[f]]))){
                                                     if (is.null(dim(object[[f]])))
                                                         p <- c(object[[f]][testdata[["id"]]])
                                                     else
                                                         p <- c(object[[f]][testdata[["id"]]])
                                                 }else{
                                                      if (!is.null(traindata)){
                                                          set.seed(trainseed)
                                                          if (f==0)
                                                              trained.model <- trainModel(model=nullobject[[1]],data=traindata)
                                                          else
                                                              trained.model <- trainModel(model=object[[f]],data=traindata)
                                                          if ("try-error" %in% class(trained.model))
                                                              stop(paste0("Failed to fit model ",f,ifelse(try(b>0,silent=TRUE),paste0(" in cross-validation step ",b,"."))))
                                                      }
                                                      else{
                                                          if (f==0)
                                                              trained.model <- nullobject[[1]]
                                                          else
                                                              trained.model <- object[[f]]
                                                      }
                                                      p <- c(do.call(predictHandlerFun, c(list(object=trained.model),args)))
                                                      if (f==0 && responseType!="binary") {## glm predicts the same value for all subjects
                                                          p <- rep(p,rep(N,NT))
                                                      }
                                                      ## predict risks not survival
                                                      if (predictHandlerFun=="predictSurvProb") p <- 1-p
                                                  }
                                                 if (!is.null(times)){
                                                     data.table(id=testdata[["id"]],model=f,risk=p,times=rep(times,rep(N,NT)))
                                                 } else {
                                                       data.table(id=testdata[["id"]],model=f,risk=p)
                                                   }
                                             }))
        if (!is.null(Weights)){
            if (Weights$method=="marginal"){
                Wt <- data.table(times=times,Wt=Weights$IPCW.times)
            } else {
                  Wt <- data.table(times=rep(times,rep(N,NT)),Wt=Weights$IPCW.times)
              }
            pred <- merge(pred,Wt,by="times")
        }
        # compute and test performance
        ## input <- list("prediction"=pred,response=response)
        input <- list(DT=merge(response,pred,by="id"),
                      N=N,
                      NT=NT,
                      NF=NF,
                      alpha=alpha,
                      test=test,
                      dolist=dolist)
        if (responseType=="competing.risks")
            input <- c(input,list(cause=cause))
        ## browser(skipCalls=1)
        if (responseType %in% c("survival","competing.risks") &&test==TRUE){
            if (censType=="rightCensored"){
                input <- c(input,list(MC=response[,getInfluenceCurve.KM(time=time,status=status)]))
            }
            else{
                input <- c(input,list(MC=matrix(0,ncol=length(response$time),nrow=length(unique(response$time)))))
            }
        }
        out <- lapply(metrics,function(m){
                          x <- do.call(paste(m,responseType,sep="."),input)
                          ## if (!is.null(x$score))
                          x$score[,model:=factor(model,levels=mlevs,mlabels)]
                          ## else
                          ## x[,model:=factor(model,levels=mlevs,mlabels)]
                          if (NROW(x$test)>0){
                              x$test[,model1:=factor(model1,levels=mlevs,mlabels)]
                              x$test[,model2:=factor(model2,levels=mlevs,mlabels)]
                          }else{
                               if (!is.null(x$test)) x$test <- NULL
                           }
                          x
                      })
        names(out) <- metrics
        out
    }
    # }}}
    # -----------------apparent nosplit performance---------------------
    # {{{
    noSplit <- computePerformance(object=object,
                                  nullobject=nullobject,
                                  testdata=data,
                                  metrics=metrics,
                                  times=times,
                                  cause=cause,
                                  test=test,
                                  alpha=alpha,
                                  dolist=dolist,
                                  NF=NF,
                                  NT=NT)
    # }}}
    # -----------------crossvalidation performance---------------------
    # {{{ 
    crossval <- NULL
    if (splitMethod$name=="BootCv"){
        if (missing(trainseeds)||is.null(trainseeds))
            trainseeds <- sample(1:1000000,size=B,replace=FALSE)
        crossval <- foreach (b=1:B) %dopar%{
                                         traindata=data[splitMethod$index[,b],,drop=FALSE]
                                         ## subset.data.table preserves order
                                         testdata <- subset(data,(match(1:N,unique(splitMethod$index[,b]),nomatch=0)==0),drop=FALSE)
                                         cb <- computePerformance(object=object,
                                                                  nullobject=nullobject,
                                                                  testdata=testdata,
                                                                  traindata=traindata,
                                                                  trainseed=trainseeds[b],
                                                                  metrics=metrics,
                                                                  times=times,
                                                                  cause=cause,
                                                                  test=test,
                                                                  alpha=alpha,
                                                                  dolist=dolist,
                                                                  NF=NF,
                                                                  NT=NT)
                                         cb
                                     }
        ## browser(skipCalls=1)
        crossvalPerf <- lapply(names(crossval[[1]]),function(m){
                                        ## if (test==TRUE){
                                        if (length(crossval[[1]][[m]]$score)>0){
                                            bootcv <- data.table::rbindlist(lapply(crossval,function(x){x[[m]]$score}))
                                            ## if (length(dolist)>0){
                                            if (length(crossval[[1]][[m]]$test)>0){
                                                if (responseType %in% c("survival","competing.risks")){
                                                    multisplit.test <- data.table::rbindlist(lapply(crossval,function(x){x[[m]]$test[,data.table(times,model1,model2,p)]}))
                                                }else{
                                                     multisplit.test <- data.table::rbindlist(lapply(crossval,function(x){x[[m]]$test[,data.table(model1,model2,p)]}))
                                                 }
                                            }else{ 
                                                 multisplit.test <- NULL
                                             }
                                        }else{
                                             multisplit.test <- NULL
                                             bootcv <- data.table::rbindlist(lapply(crossval,function(x){x[[m]]}))
                                         }
                                        if (responseType %in% c("survival","competing.risks")){
                                            bootcv <- bootcv[,mean(eval(as.name(m))),by=list(model,times)]
                                            data.table::setnames(bootcv,c("model","times",m))
                                        } else{
                                              bootcv <- bootcv[,mean(eval(as.name(m))),by=list(model)]
                                              data.table::setnames(bootcv,c("model",m))
                                          }
                                        out <- list(bootcv)
                                        names(out) <- "score"
                                        ## names(out) <- paste0("Cross-validation (average of ",B," steps)")
                                        if (!is.null(multisplit.test)){
                                            ms <- list(multisplit.test)
                                            names(ms) <- "tests"
                                            ## names(ms) <- paste0("Multisplit test (",B," splits)")
                                            out <- c(out,ms)
                                        }
                                        out
                                    })
        names(crossvalPerf) <- names(crossval[[1]])
        ## Brier <- data.table::rbindlist(lapply(crossval,function(x)x[["Brier"]]))
        ## Brier[,list(looboot=mean(ipcwResiduals)),by=list(model,times,id)]
        ## bootcv=Brier[,list(bootcv=mean(ipcwResiduals)),by=list(model,times,b)]
    }
    # }}}
    #------------------output-----------------------------------
    if (is.null(crossval))
        output <- noSplit
    else{
        ## output <- list(noSplitPerf=noSplit,crossValPerf=crossvalPerf)
        output <- crossvalPerf
    }
    # -----------------specific task handlers-------------------
    # {{{ 
    output <- c(output,list(responseType=responseType,
                            predictHandlerFun=predictHandlerFun,
                            censType=censType,
                            censoringHandling=censMethod,
                            splitMethod=splitMethod,metrics=metrics))
    for (m in metrics){
        output[[m]]$metric <- m
        class(output[[m]]) <- "score"
    }
    class(output) <- "Score"
    output
}

#' @export
print.Score <- function(x,...){
    B <- x$splitMethod$B
    for (m in x$metrics){
        cat(paste0("\nMetric ",m,":\n"))
        print(x[[m]],B)
    }
}
#' @export
print.score <- function(x,B=0,...){
    if (B>0){
        print(x$score)
        cat(paste0("\nCross-validation (average of ",B," steps)\n\n"))
        if (!is.null(x$tests)){
            print(x$tests)
            cat(paste0("\nMultisplit test (",B," splits)\n\n"))
        }
    }else{
         print(x$score)
         if (!is.null(x$tests))
             print(x$tests)
     }
}

Brier.binary <- function(DT,test=TRUE,alpha,N,NT,NF,dolist){
    residuals=risk=model=ReSpOnSe=lower.Brier=upper.Brier=se.Brier=NULL
    DT[,residuals:=(ReSpOnSe-risk)^2,by=model]
    score <- DT[,list(Brier=mean(residuals)),by=model]
    if (test==TRUE){
        data.table::setorder(DT,model,ReSpOnSe)
        score <- DT[,data.table(Brier=sum(residuals)/N,se.Brier=sd(residuals)/sqrt(N)),by=list(model)]
        score[,lower.Brier:=Brier-qnorm(1-alpha/2)*se.Brier]
        score[,upper.Brier:=Brier + qnorm(1-alpha/2)*se.Brier]
        data.table::setkey(score,model)
        data.table::setkey(DT,model)
        DT <- DT[score]
        test.Brier <- DT[,getComparisons(data.table(x=Brier,IC=residuals,model=model),NF=NF,N=N,alpha=alpha,dolist=dolist)]
        Brier <- list(score=score,test=test.Brier)
    }else{
         Brier <- list(score=DT[,list(Brier=mean(residuals)),by=list(model)])
     }
    Brier
}
auRoc.numeric <- function(X,D,breaks){
    if (is.null(breaks)) breaks <- rev(sort(unique(X))) ## need to reverse when high X is concordant with {response=1}  
    TPR <- c(prodlim::sindex(jump.times=X[D==1],eval.times=breaks,comp="greater",strict=FALSE)/sum(D==1),0)
    FPR <- c(prodlim::sindex(jump.times=X[D==0],eval.times=breaks,comp="greater",strict=FALSE)/sum(D==0),0)
    0.5 * sum(diff(c(0,FPR,1)) * (c(TPR,1) + c(0,TPR)))
}
auRoc.factor <- function(X,D){
    TPR <- (sum(D==1)-table(X[D==1]))/sum(D==1)
    FPR <- table(X[D==0])/sum(D==0)
    0.5 * sum(diff(c(0,FPR,0,1)) * (c(TPR,0,1) + c(0,TPR,0)))
}
AUC.binary <- function(DT,breaks=NULL,test,alpha,N,NT,NF,dolist){
    model=risk=ReSpOnSe=NULL
    data.table::setorder(DT,model,risk)
    if (is.factor(DT[["risk"]])){
        score <- DT[,list(AUC=auRoc.factor(risk,ReSpOnSe)),by=list(model)]
    }
    else{
        score <- DT[,list(AUC=auRoc.numeric(risk,ReSpOnSe,breaks=NULL)),by=list(model)]
    }
    AUC <- list(score=score)
    AUC
}
Brier.survival <- function(DT,MC,test,alpha,N,NT,NF,dolist){
    Yt=time=times=Residuals=risk=ipcwResiduals=WTi=Wt=status=setorder=model=IC.Brier=data.table=sd=lower.Brier=qnorm=se.Brier=upper.Brier=NULL
    ## compute 0/1 outcome:
    DT[,Yt:=1*(time<=times)]
    ## compute residuals
    DT[,Residuals:=(Yt-risk)^2]
    ## apply weights 
    DT[,ipcwResiduals:=Residuals/WTi]
    DT[Yt==0,ipcwResiduals:=Residuals/Wt]
    ## deal with censored observations before t
    DT[Yt==1 & status==0,ipcwResiduals:=0]
    ## DT[,c("Yt","time","WTi","Wt","status","Residuals"):=NULL]
    if (test==TRUE){
        data.table::setorder(DT,model,times,time,-status)
        DT[,IC.Brier:=getInfluenceCurve.Brier(t=times[1],time=time,Yt=Yt,ipcwResiduals=ipcwResiduals,MC=MC),by=list(model,times)]
        score <- DT[,data.table(Brier=sum(ipcwResiduals)/N,se.Brier=sd(IC.Brier)/sqrt(N)),by=list(model,times)]
        score[,lower.Brier:=Brier-qnorm(1-alpha/2)*se.Brier]
        score[,upper.Brier:=Brier + qnorm(1-alpha/2)*se.Brier]
        data.table::setkey(score,model,times)
        data.table::setkey(DT,model,times)
        DT <- DT[score]
        test.Brier <- DT[,getComparisons(data.table(x=Brier,IC=IC.Brier,model=model),NF=NF,N=N,alpha=alpha,dolist=dolist),by=list(times)]
        Brier <- list(score=score,test=test.Brier)
    }else{
         Brier <- list(score=DT[,data.table(Brier=sum(ipcwResiduals)/N),by=list(model,times)])
     }
    Brier
}
Brier.competing.risks <- function(DT,MC,test,alpha,N,NT,NF,dolist,cause){
    Yt=time=times=event=Residuals=risk=ipcwResiduals=WTi=Wt=status=setorder=model=IC.Brier=data.table=sd=lower.Brier=qnorm=se.Brier=upper.Brier=NULL
    ## compute 0/1 outcome:
    DT[,Yt:=1*(time<=times & event==cause)]
    ## compute residuals
    DT[,Residuals:=(Yt-risk)^2]
    ## apply weights 
    DT[,ipcwResiduals:=Residuals/WTi]
    DT[time>times,ipcwResiduals:=Residuals/Wt]
    ## deal with censored observations before t
    DT[time<=times & status==0,ipcwResiduals:=0]
    ## DT[,c("Yt","time","WTi","Wt","status","Residuals"):=NULL]
    if (test==TRUE){
        data.table::setorder(DT,model,times,time,-status)
        DT[,IC.Brier:=getInfluenceCurve.Brier(t=times[1],time=time,Yt=Yt,ipcwResiduals=ipcwResiduals,MC=MC),by=list(model,times)]
        score <- DT[,data.table(Brier=sum(ipcwResiduals)/N,se.Brier=sd(IC.Brier)/sqrt(N)),by=list(model,times)]
        score[,lower.Brier:=Brier-qnorm(1-alpha/2)*se.Brier]
        score[,upper.Brier:=Brier + qnorm(1-alpha/2)*se.Brier]
        data.table::setkey(score,model,times)
        data.table::setkey(DT,model,times)
        DT <- DT[score]
        test.Brier <- DT[,getComparisons(data.table(x=Brier,IC=IC.Brier,model=model),NF=NF,N=N,alpha=alpha,dolist=dolist),by=list(times)]
        Brier <- list(score=score,test=test.Brier)
    }else{
         Brier <- list(score=DT[,data.table(Brier=sum(ipcwResiduals)/N),by=list(model,times)])
     }
    Brier
}

AireTrap <- function(FP,TP,N){
    N <- length(FP)
    sum((FP-c(0,FP[-N]))*((c(0,TP[-N])+TP)/2))
}

AUC.survival <- function(DT,MC,test,alpha,N,NT,NF,dolist){
    model=times=risk=Cases=time=status=Controls1=TPt=FPt=WTi=Wt=ipcwControls1=ipcwCases=IC.AUC=lower.AUC=se.AUC=upper.AUC=NULL
    cause <- 1
    ## assign Weights before ordering
    DT[,ipcwControls1:=1/(Wt*N)]
    ## DT[,ipcwControls2:=1/(WTi*N)]
    DT[,ipcwCases:=1/(WTi*N)]
    ## DT[,ipcwControls2:=1/(WTi*N)]
    ## order data
    data.table::setorder(DT,model,times,-risk)
    ## identify cases and controls
    DT[,Cases:=(time <= times &  status==cause)]
    DT[,Controls1:=(time > times)] 
    ## DT[,Controls2:=(time < times &  status!=cause & status !=0)]
    ## prepare Weights
    DT[Cases==0,ipcwCases:=0]
    DT[Controls1==0,ipcwControls1:=0]
    ## DT[Controls2==0,ipcwControls2:=0]
    ## compute denominator
    ## ROC <- DT[,list(TPt=c(0,cumsum(ipcwCases)),FPt=c(0,cumsum(ipcwControls1)+cumsum(ipcwControls2))),by=list(model,times)]
    DT[,TPt:=cumsum(ipcwCases)/sum(ipcwCases),by=list(model,times)]
    ## DT[,FPt:=(cumsum(ipcwControls1)+cumsum(ipcwControls2))/(sum(ipcwControls2)+sum(ipcwControls1)),by=list(model,times)]
    DT[,FPt:=(cumsum(ipcwControls1))/(sum(ipcwControls1)),by=list(model,times)]
    nodups <- DT[,c(!duplicated(risk)[-1],TRUE),by=list(model,times)]$V1
    ## auc <- DT[nodups,list(AUC=sum((FPt-c(0,FPt[-N]))*((c(0,TPt[-N])+TPt)/2))),by=list(model,times)]
    score <- DT[nodups,list(AUC=AireTrap(FPt,TPt)),by=list(model,times)]
    data.table::setkey(score,model,times)
    if (test==TRUE){
        ## compute influence function
        data.table::setorder(DT,model,times,time,-status)
        DT[,IC.AUC:=getInfluenceCurve.AUC.survival(t=times[1],n=N,time=time,risk=risk,Cases=Cases,Controls1=Controls1,ipcwControls1=ipcwControls1,ipcwCases=ipcwCases,MC=MC), by=list(model,times)]
        se.score <- DT[,list(se.AUC=sd(IC.AUC)/sqrt(N)),by=list(model,times)]
        data.table::setkey(se.score,model,times)
        score <- score[se.score]
        score[,lower.AUC:=AUC-qnorm(1-alpha/2)*se.AUC]
        score[,upper.AUC:=AUC+qnorm(1-alpha/2)*se.AUC]
        data.table::setkey(DT,model,times)
        DT <- DT[score]
        test.AUC <- DT[,getComparisons(data.table(x=AUC,IC=IC.AUC,model=model),NF=NF,N=N,alpha=alpha,dolist=dolist),by=list(times)]
        AUC <- list(score=score,test=test.AUC)
    }else{
         AUC <- list(score=score)
     }
    AUC
}


AUC.competing.risks <- function(DT,MC,test,alpha,N,NT,NF,dolist,cause){
    model=times=risk=Cases=time=status=event=Controls1=Controls2=TPt=FPt=WTi=Wt=ipcwControls1=ipcwControls2=ipcwCases=IC.AUC=lower.AUC=se.AUC=upper.AUC=NULL
    ## assign Weights before ordering
    DT[,ipcwControls1:=1/(Wt*N)]
    DT[,ipcwControls2:=1/(WTi*N)]
    DT[,ipcwCases:=1/(WTi*N)]
    DT[,ipcwControls2:=1/(WTi*N)]
    ## order data
    data.table::setorder(DT,model,times,-risk)
    ## identify cases and controls
    DT[,Cases:=(time <=times &  event==cause)]
    DT[,Controls1:=(time > times)] 
    DT[,Controls2:=(time <=times &  event!=cause & status !=0)]
    ## prepare Weights
    DT[Cases==0,ipcwCases:=0]
    DT[Controls1==0,ipcwControls1:=0]
    DT[Controls2==0,ipcwControls2:=0]
    ## compute denominator
    ## ROC <- DT[,list(TPt=c(0,cumsum(ipcwCases)),FPt=c(0,cumsum(ipcwControls1)+cumsum(ipcwControls2))),by=list(model,times)]
    DT[,TPt:=cumsum(ipcwCases)/sum(ipcwCases),by=list(model,times)]
    DT[,FPt:=(cumsum(ipcwControls1)+cumsum(ipcwControls2))/(sum(ipcwControls2)+sum(ipcwControls1)),by=list(model,times)]
    nodups <- DT[,c(!duplicated(risk)[-1],TRUE),by=list(model,times)]$V1
    ## auc <- DT[nodups,list(AUC=sum((FPt-c(0,FPt[-N]))*((c(0,TPt[-N])+TPt)/2))),by=list(model,times)]
    score <- DT[nodups,list(AUC=AireTrap(FPt,TPt)),by=list(model,times)]
    data.table::setkey(score,model,times)
    if (test==TRUE){
        ## compute influence function
        data.table::setorder(DT,model,times,time,-status)
        DT[,IC.AUC:=getInfluenceCurve.AUC.competing.risks(t=times[1],n=N,time=time,risk=risk,ipcwControls1=ipcwControls1,ipcwControls2=ipcwControls2,ipcwCases=ipcwCases,Cases=Cases,Controls1=Controls1,Controls2=Controls2,MC=MC), by=list(model,times)]
        se.score <- DT[,list(se.AUC=sd(IC.AUC)/sqrt(N)),by=list(model,times)]
        data.table::setkey(se.score,model,times)
        score <- score[se.score]
        score[,lower.AUC:=AUC-qnorm(1-alpha/2)*se.AUC]
        score[,upper.AUC:=AUC+qnorm(1-alpha/2)*se.AUC]
        data.table::setkey(DT,model,times)
        DT <- DT[score]
        test.AUC <- DT[,getComparisons(data.table(x=AUC,IC=IC.AUC,model=model),NF=NF,N=N,alpha=alpha,dolist=dolist),by=list(times)]
        AUC <- list(score=score,test=test.AUC)
    }else{
         AUC <- list(score=score)
     }
    AUC
}


