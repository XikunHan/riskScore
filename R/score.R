##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param object 
##' @param ... 
##' @return 
##' @author Thomas Alexander Gerds
score <- function(object,...){
  UseMethod("score",object=object)
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param object 
##' @param formula 
##' @param data 
##' @param event 
##' @param metrics 
##' @param times 
##' @param maxtime 
##' @param landmarks 
##' @param useEventTimes 
##' @param nullModel 
##' @param censMethod 
##' @param censModel 
##' @param splitMethod 
##' @param B 
##' @param M 
##' @param verbose 
##' @param ... 
##' @return 
##' @author Thomas Alexander Gerds
score.list <- function(object,
                       formula,
                       data,
                       event,
                       metrics=c("auc","brier"),
                       tests=c("t.test"),
                       plots=c("roc","calibration"),
                       times,
                       maxtime,
                       landmarks,
                       useEventTimes=FALSE,
                       nullModel=TRUE,
                       censMethod="ipcw",
                       censModel="cox",
                       splitMethod,
                       B,
                       M,
                       verbose=0,
                       trainseeds,
                       ...){
    theCall <- match.call(expand=TRUE)
    # -----------------parse arguments and prepare data---------
    # {{{ Response

    formula <- getFormula(formula,object,verbose=verbose)
    data <- data.table(data)
    ## getData <- function(data,formula,object,verbose,...){
    ## if (missing(data)){
    ## data <- eval(object[[1]]$call$data)
    ## if (match("data.frame",class(data),nomatch=0)==0)
    ## stop("Argument data is missing.")
    ## else
    ## if (verbose)
    ## warning("Argument data is missing. I use the data from the call to the first model instead.")
    ## }
    ## data
    ## }
    ## data <- getData(data=data,
    ## formula=formula,
    ## object=object,
    ## verbose=verbose)
    responseFormula <- update(formula,~1)
    if (missing(event)) event <- NULL
    responsevars <- all.vars(responseFormula)
    response <- getResponse(formula=responseFormula,
                            data=data,
                            event=event,
                            vars=responsevars)
    response.dim <- NCOL(response)
    responseType <- attr(response,"model")
    ## make sure response variable names are time, event, status
    ## and that response is in the first column(s) 
    data[,eval(responsevars):=NULL]
    data <- cbind(response,data)
    if (responseType=="survival")
        formula <- update(formula,"Hist(time,status)~.")
    if (responseType=="competing.risks")
        formula <- update(formula,"Hist(time,event)~.")
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
    allmethods <- methods(predictHandlerFun)
    ## wantedMethods <- lapply(object,function(o){
    ## candidateMethods <- paste(predictHandlerFun,class(o),sep=".")
    ## if (all(match(candidateMethods,allmethods,nomatch=0)==0))
    ## stop(paste("Could not find ",predictHandlerFun," method for ",paste(class(o),collapse=" ,"),sep=""))
    ## })
    # checking the models for compatibility with resampling
    if (do.resample){
        object <- lapply(1:length(object),function(f){
                                     fit <- object[[f]]
                                     if(is.null(fit$call))
                                         stop(paste("model",names(object)[f],"does not have a call argument."))
                                     else fit$call$data <- NULL
                                     fit
                                 })
    }
    # add null model and find names for the object
    if (nullModel==TRUE)
        object <- c(getNullModel(formula=formula,
                                 data=data,
                                 responseType=responseType),object)
    if (is.null(names(object))){
        names(object) <- sapply(object,function(o)class(o)[1])}
    else {
        names(object)[(names(object)=="")] <- sapply(object[(names(object)=="")],function(o)class(o)[1])
    }
    names(object) <- make.names(names(object),unique=TRUE)
    NF <- length(object)
    ## if (any(whichNum <- sapply(sapply(object,class),function(cc)sum(match(cc,c("numeric","matrix"),nomatch=0))))>0)
    ## warning("Numeric or matrix prediction. This works only if data and object have the same order.",call.=FALSE)
    # }}}
    # {{{ add id *before* ordering the data. data have to be ordered for pec::ipcw
    data[,id:=1:N]
    # }}}
    # {{{ Evaluation landmarks and horizons (times)
    if (responseType %in% c("survival","competing.risks")){
        ## order according to convention, in case of a tie,
        ## events come before right censored
        setorder(data,time,-status)
        if (missing(maxtime) || is.null(maxtime)){
            eventTimes <- unique(data[,time])
            maxtime <- eventTimes[length(eventTimes)]
        }
        if (missing(landmarks)){
            start <- 0
            if (missing(times)){
                if (useEventTimes==TRUE)
                    times <- unique(c(start,eventTimes))
                else
                    times <- seq(start,maxtime,(maxtime - start)/100)
            }
            else{
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
        if (!missing(times)) warning("Function 'score': Response type is not time-to-event: argument 'times' will be ignored.",call.=FALSE)
        times <- NULL
    }
    # }}}
    # ----------------------------find metrics and plots ----------------------
    # {{{
    
    ## Metrics <- lapply(metrics,grep,c("AUC","Brier"),ignore.case=TRUE,value=TRUE)
    Metrics <- metrics
    Plots <- lapply(plots,grep,c("Roc","Cal"),ignore.case=TRUE,value=TRUE)
    # }}}
    # -----------------specific task handlers-------------------
    # {{{ 

    if (verbose==1){
        cat(paste("\nresponseType:",responseType))
        cat(paste("\npredictHandlerFun:",predictHandlerFun))
        cat(paste("\ncensType:",censType))
        cat(paste("\ncensoringHandling:",censMethod))
        cat(paste("\nsplitMethod:",splitMethod$name))
        cat("\n")
    }
    # }}}
    # -----------------IPCW outside loop -----------------------
    # {{{ 
    if (responseType %in% c("survival","competing.risks")){
        if (censType=="rightCensored"){
            if ("outside" %in% censMethod){
                if ("ipcw" %in% censMethod){
                    Weights <- getCensoringWeights(formula=formula,data=testdata,times=times,censModel=censModel,responseType=responseType)
                }
                else{
                    censMethod <- "jackknife.pseudo.values"
                    pseudoResponse <- getPseudoValues(formula=formula,data=testdata,responseType=responseType,times=times,event=event)
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
    computePerformance <- function(traindata,
                                   testdata,
                                   trainseed,
                                   metrics,
                                   tests){
        NF <- length(object)
        # split into response and predictors
        response <- testdata[,1:response.dim,with=FALSE]
        response[,id:=testdata[["id"]]]
        X <- testdata[,-c(1:response.dim),with=FALSE]
        # {{{ -----------------IPCW inner loop-----------------------
        if ("inside" %in% censMethod){
            if (censType=="rightCensored"){
                if ("ipcw" %in% censMethod){
                    Weights <- getCensoringWeights(formula=formula,
                                                   data=testdata,
                                                   times=times,
                                                   censModel=censModel,
                                                   responseType=responseType)
                    response[,WTi:=Weights$IPCW.subjectTimes]
                } else{
                      censMethod <- "jackknife.pseudo.values"
                      pseudoResponse <- getPseudoValues(formula=formula,data=testdata,responseType=responseType,times=times,event=event)
                  }
            } else{
                  if (censType=="uncensored"){
                      Weights <- list(IPCW.times=rep(1,NT),
                                      IPCW.subjectTimes=matrix(1,ncol=NT,nrow=NROW(testdata)))
                      response[,WTi:=1]
                  } else{
                        stop("Cannot handle this type of censoring.")
                    }
              }
        }
        # }}}
        # extract predictions as data.table
        args <- switch(predictHandlerFun,"predictStatusProb"={list(newdata=X)},
                       "predictSurvProb"={list(newdata=X,times=times)},
                       "predictEventProb"={list(newdata=X,times=times,event=event)},
                       stop("Unknown predictHandler."))
        ## if (!is.null(traindata)){
        ## set.seed(trainseed)
        ## trained.nullmodel <- trainModel(model=object[[1]],data=traindata)
        ## } else{
        ## trained.nullmodel <- object[[1]]
        ## }
        ## predNull <- data.table(do.call(predictHandlerFun, c(list(object=trained.nullmodel),args)))
        ## print(NROW(testdata))
        ## print(NROW(traindata))
        pred <- rbindlist(lapply(2:NF,function(f){
                                     if (any(c("integer","factor","numeric","matrix") %in% class(object[[f]]))){
                                         if (is.null(dim(object[[f]])))
                                             p <- c(object[[f]][testdata[,id]])
                                         else
                                             p <- c(object[[f]][testdata[,id],])
                                     }else{
                                          if (!is.null(traindata)){
                                              set.seed(trainseed)
                                              trained.model <- trainModel(model=object[[f]],data=traindata)
                                          }
                                          else{
                                              trained.model <- object[[f]]
                                          }
                                          p <- c(do.call(predictHandlerFun, c(list(object=trained.model),args)))
                                          ## predict risks not survival
                                          if (predictHandlerFun=="predictSurvProb") p <- 1-p
                                      }
                                     if (!is.null(times)){
                                         if (!is.null(Weights)){
                                             if (Weights$method=="marginal")
                                                 Wt <- rep(Weights$IPCW.times,rep(NROW(testdata),NT))
                                             else
                                                 Wt <- c(Weights$IPCW.times)
                                             data.table(id=testdata[["id"]],
                                                        model=f,
                                                        pred=p,
                                                        times=rep(times,rep(NROW(testdata),NT)),
                                                        Wt=Wt)
                                         } else{
                                               data.table(id=testdata[["id"]],model=f,pred=p,times=rep(times,rep(N,NT)))
                                           }
                                     }else{
                                          data.table(id=testdata[["id"]],model=f,pred=p)
                                      }
                                 }))
        # compute and test performance
        ## input <- list("prediction"=pred,response=response)
        input <- merge(response,pred,by="id")
        ## if (!is.null(Weights)) input <- c(input,list(Weights=Weights))
        ## if (!is.null(Wt)) input <- c(input,list(Weights=Wt))
        ## if (!is.null(times)) input <- c(input,list(times=times))        
        setattr(input,"cause",1)
        if (!is.null(times)) setattr(input,"NT",NT)
        setattr(input,"NF",NF)
        setattr(input,"N",NROW(testdata))
        out <- lapply(metrics,function(m){do.call(paste(m,responseType,sep="."),list(input))})
        names(out) <- metrics
        out
    }
    # }}}
    # -----------------apparent nosplit performance---------------------
    # {{{
    noSplit <- computePerformance(traindata=NULL,
                                  testdata=data,
                                  trainseed=NULL,
                                  metrics=Metrics)
    crossval <- NULL
    if (splitMethod$name=="BootCv"){
        if (missing(trainseeds)||is.null(trainseeds))
            trainseeds <- sample(1:1000000,size=B,replace=FALSE)
        if (require(foreach)){
            crossval <- foreach (b=1:B) %dopar% {
                                             ## print(b)
                                             train=data[splitMethod$index[,b],,drop=FALSE]
                                             ## subset.data.table preserves order
                                             test <- subset(data,(match(1:N,unique(splitMethod$index[,b]),nomatch=0)==0),drop=FALSE)
                                             cb <- computePerformance(b,
                                                                      traindata=train,
                                                                      trainseed=trainseeds[b],
                                                                      testdata=test,
                                                                      metrics=Metrics)
                                             ## browser(skipCalls=1)
                                             cb$Brier[,b:=b]
                                             cb
                                         }
        }
        ## tmp <- lapply(Metrics,function(m)rbindlist(lapply(crossval,function(x)x[[m]])))
        Brier <- rbindlist(lapply(crossval,function(x)x[["Brier"]]))
        Brier[,list(looboot=mean(Rtw)),by=list(model,times,id)]
        bootcv=Brier[,list(bootcv=mean(Rtw)),by=list(model,times,b)]
    }
    # }}}
    #------------------output-----------------------------------
    ## browser(skipCalls=1)
    output <- list(noSplitPerf=noSplit,crossValPerf=crossval)
}

Brier.binary <- function(prediction,response){
    prediction[,residuals:=(response-p1)^2,by=model]
    bs <- prediction[,list(Brier=mean(residuals)),by=model]
    bs
}

## AUC.binary <- function(prediction,response){
## auc <- prediction[model!=1,list(AUC=auc.binary(p1,response)),by=model]
## auc    
## }

auRoc.numeric <- function(X,D,breaks){
    if (is.null(breaks)) breaks <- sort(unique(X))
    TPR <- prodlim::sindex(jump.times=X[D==1],eval.times=breaks,comp="greater")/sum(D==1)
    FPR <- prodlim::sindex(jump.times=X[D==0],eval.times=breaks,comp="greater")/sum(D==0)
    ## FIXME: maybe add 0 at end of TPR and FPR at earlier stage?
    ## 0.5 * sum(diff(c(0,FPR,1)) * (c(TPR,1) + c(0,TPR)))
    0.5 * sum(diff(c(0,FPR,0,1)) * (c(TPR,0,1) + c(0,TPR,0)))
}
auRoc.factor <- function(X,D){
    TPR <- (sum(D==1)-table(X[D==1]))/sum(D==1)
    FPR <- table(X[D==0])/sum(D==0)
    ## FIXME: maybe add 0 at end of TPR and FPR at earlier stage?
    ## 0.5 * sum(diff(c(0,FPR,1)) * (c(TPR,1) + c(0,TPR)))
    0.5 * sum(diff(c(0,FPR,0,1)) * (c(TPR,0,1) + c(0,TPR,0)))
}

AUC.binary <- function(object,breaks=NULL){
    if (is.factor(object[["pred"]])){
        object[,auRoc.factor(pred,Y),by=list(model)]
    }
    else{
        object[,auRoc.numeric(pred,Y,breaks=NULL),by=list(model)]
    }
}
Brier.survival <- function(object,se=TRUE,residuals=FALSE){
    n <- object[model==2&times==min(times),length(time)]
    ## compute 0/1 outcome:
    ## object[,Yt:=1*(time>=times)]
    object[,Yt:=1*(time<=times)]
    ## compute residuals
    object[,Rt:=(Yt-pred)^2]
    ## apply weights 
    ## prediction[Yt==0,Rtw:=response[["WTi"]]*Rt]
    object[,Rtw:=Rt/WTi]
    object[Yt==0,Rtw:=Rt/Wt]
    ## deal with censored observations before t
    object[Yt==1 & status==0,Rtw:=0]
    ##
    if (se==TRUE){
        setorder(object,model,times,time,-status)
        MatInt0TcidhatMCksurEff  <- object[model==2&times==min(times),getInfluenceCurve.KM(time=time,status=status)]
        getInfluenceCurve.BS(object,MatInt0TcidhatMCksurEff)
    }
    object[,c("Yt","time","WTi","Wt","status","Rt"):=NULL]
    ## tests <- lapply(unique(prediction$times),function(t){
    ## prediction[times==t,pairwise.t.test(x=Rtw,g=model)]})
    if (residuals==TRUE)
        object
    else
        object[,data.table(Brier=sum(Rtw)/n,se=sd(IC)/sqrt(n)),by=list(model,times)]
}

Brier1.survival <- function(prediction,response,times,Weights){
    ## we could merge prediction and response by id
    ## for now we use that NROW(prediction) is a multiple 
    ## of NROW(response)    
    ## 
    n <- NROW(response)
    ## compute 0/1 outcome 
    prediction[,Yt:=1*(response[["time"]]>=times)]
    ## compute residuals
    prediction[,Rt:=(Yt-pred)^2]
    ## apply weights 
    ## prediction[Yt==0,Rtw:=response[["WTi"]]*Rt]
    prediction[,Rtw:=Rt/response[["WTi"]]]
    prediction[Yt==1,Rtw:=Rt/Wt]
    ## get rid of censored observations before time. 
    ## make this more memory efficient see ~/Howto/data.table.org
    prediction <- prediction[Yt==1 | response[["status"]]!=0]
    ## browser(skipCalls=1)
    ## tests <- lapply(unique(prediction$times),function(t){
    ## prediction[times==t,pairwise.t.test(x=Rtw,g=model)]})
    prediction[,data.table(Brier=sum(Rtw)/n),by=list(model,times)]
}
Brier2.survival <- function(prediction,response,times,Weights){
    NT <- length(times)
    N <- NROW(response)
    NF <- max(prediction$model)
    ## browser(skipCalls=1)
    bs <- rbindlist(lapply(2:NF,function(f){
                               data.table(model=f,times=times,Brier=.C("pec",
                                                                  pec=double(NT),
                                                                  as.double(response[["time"]]),
                                                                  as.double(response[["status"]]),
                                                                  as.double(times),
                                                                  as.double(unlist(c(prediction[model==f,"pred",with=FALSE]))),
                                                                  as.double(Weights),
                                                                  as.double(response[["WTi"]]),
                                                                  as.integer(N),
                                                                  as.integer(NT),
                                                                  as.integer(0),
                                                                  ## as.integer(is.null(dim(prediction))),
                                                                  as.integer(0),
                                                                  NAOK=TRUE,
                                                                  PACKAGE="pec")$pec)
                           }))
    bs
}

# fonction: area under the curve by trapezoidal rule
AireTrap<-function(FP,TP){
    N <- length(FP)
    sum((FP-c(0,FP[-N]))*((c(0,TP[-N])+TP)/2))
}

AUC.survival <- function(Matdata){
    cause <- attr(Matdata,"cause")
    NT <- attr(Matdata,"NT")
    NF <- attr(Matdata,"NF")
    N <- attr(Matdata,"N")
    cause <- 1
    ## assign Weights before ordering
    Matdata[,ipcwControls1:=1/(Wt*N)]
    Matdata[,ipcwControls2:=1/(WTi*N)]
    Matdata[,ipcwCases:=1/(WTi*N)]
    Matdata[,ipcwControls2:=1/(WTi*N)]
    ## order data
    setorder(Matdata,model,times,pred)
    ## identify cases and controls
    Matdata[,Cases:=(time < times &  status==cause)]
    Matdata[,Controls1:=(time > times)] 
    Matdata[,Controls2:=(time < times &  status!=cause & status !=0)]
    ## prepare Weights
    Matdata[Cases==0,ipcwCases:=0]
    Matdata[Controls1==0,ipcwControls1:=0]
    Matdata[Controls2==0,ipcwControls2:=0]
    ## compute denominator
    ## ROC <- Matdata[,list(TPt=c(0,cumsum(ipcwCases)),FPt=c(0,cumsum(ipcwControls1)+cumsum(ipcwControls2))),by=list(model,times)]
    Matdata[,TPt:=cumsum(ipcwCases)/sum(ipcwCases),by=list(model,times)]
    Matdata[,FPt:=(cumsum(ipcwControls1)+cumsum(ipcwControls2))/(sum(ipcwControls2)+sum(ipcwControls1)),by=list(model,times)]
    nodups <- Matdata[,c(!duplicated(pred)[-1],TRUE),by=list(model,times)]$V1
    auc <- Matdata[nodups,list(AUC=AireTrap(FPt,TPt)),by=list(model,times)]
    setkey(auc,model,times)
    ## compute se.fit
    setorder(Matdata,model,times,time,-status)
    MatInt0TcidhatMCksurEff  <- Matdata[model==2&times==min(times),getInfluenceCurve.KM(time=time,status=status)]
    seAuc <- Matdata[,list(se=getInfluenceCurve.AUC.survival(data.table(time,status,pred,times,ipcwControls1,ipcwCases,Cases,Controls1), MatInt0TcidhatMCksurEff=MatInt0TcidhatMCksurEff)),by=list(model,times)]
    setkey(seAuc,model,times)
    auc[seAuc]
}



