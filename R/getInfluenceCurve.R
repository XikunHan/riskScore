## CAUTION : data should be oredered by time
getInfluenceCurve.AUC.survival <- function(object,
                                           MatInt0TcidhatMCksurEff){  
    Matdata <- object
    t <- Matdata[,times[1]]
    n <- Matdata[,length(time)]
    Cases <- Matdata[["Cases"]]
    Controls1 <- Matdata[["Controls1"]]
    F01t <- Matdata[,sum(ipcwCases)]
    St <- Matdata[,sum(ipcwControls1)]
    # compute frequencies of cases and controls to define 
    #the size of the matrix  Mathtij1 
    nbCases <- sum(Cases)
    nbControls1 <- sum(Controls1)
    mcase <- matrix(Matdata[Cases,pred],nrow=nbCases,ncol=nbControls1)
    mcontrol <- matrix(Matdata[Controls1,pred],nrow=nbCases,ncol=nbControls1,byrow=TRUE)
    wcase <- matrix(Matdata[Cases,ipcwCases],nrow=nbCases,ncol=nbControls1)
    wcontrol <- matrix(Matdata[Controls1,ipcwControls1],nrow=nbCases,ncol=nbControls1,byrow=TRUE)
    Mathtij1 <- (1*(mcase>mcontrol)+.5*(mcase==mcontrol))*wcase*wcontrol*n*n
    ht <- (sum(Mathtij1))/(n*n) 
    vectdit <- Matdata[["Cases"]]*Matdata[["ipcwCases"]]*n
    ## print(c(mean(vectdit),F01t))
    # Compute the vecor of all sum_{i=1}^n of {\hat{h}_{tij}}_1 for all j
    colSumsMathtij1 <- rep(0,n) # initialise at 0
    colSumsMathtij1[Cases] <- rowSums(Mathtij1) # when i is a case,  then we sum the column of  Mathtij1  
    # Compute the vecor of all sum_{j=1}^n of {\hat{h}_{tij}}_1 for all i  
    rowSumsMathtij1 <- rep(0,n) # initialize at 0
    rowSumsMathtij1[Controls1] <- colSums(Mathtij1)# when  j is a control 1, then we sum the row of  Mathtij1
    hathtstar <- (sum(Mathtij1))/(n*n)  
    # compute the vector of \frac{1_{\tilde{time}_i>=t}}{ \hat{S}_{\tilde{time}}(t)}
    vectTisupt <- n*Controls1/sum(Controls1)
    T1 <- colSums(crossprod(Mathtij1,1+MatInt0TcidhatMCksurEff[Cases,]))/n
    T3 <- colSums(hathtstar*(vectTisupt + (vectdit*(1+MatInt0TcidhatMCksurEff)-F01t)/F01t))
    ## sumijakfixe <- function(k){
    ## term1 <- Mathtij1*(1+MatInt0TcidhatMCksurEff[Cases,k]) 
    ## term2 <- vectdit*(1+MatInt0TcidhatMCksurEff[,k])
    ## term3 <- (hathtstar)*(vectTisupt + (1/F01t)*(term2-F01t))
    ## sum(term1)/n - sum(term3) 
    ## }
    ## Lessumijakfixe <- numeric(n)
    ## for (i in 1:n) Lessumijakfixe[i] <- sumijakfixe(i)/(F01t*St)  
    ## print(round(Lessumijakfixe,4))
    ## print(round((T1-T3)/(F01t*St),4))
    Lessumijakfixe <- (T1-T3)/(F01t*St)
    Lessumikajfixe <- (rowSumsMathtij1 - n*hathtstar)/(F01t*St)
    Lessumjkaifixe <-  (colSumsMathtij1 - n*hathtstar*(vectTisupt+(1/F01t)*(vectdit-F01t)))/(F01t*St)
    # We compute the iid representation of the AUC estimator
    hatIFstar <-  (Lessumijakfixe + Lessumikajfixe +  Lessumjkaifixe)/(n)
    # }}}
    # we compute the standard error of the AUC estimator
    seAUC <- sd(hatIFstar)/sqrt(n)
    ## list(iidrepresentationAUCstar=hatIFstar,
    ## seAUC=seAUC)
    seAUC
}

getInfluenceCurve.BS <- function(object,MatInt0TcidhatMCksurEff){ 
    t <- object[,times[1]]
    n <- object[,length(time)]
    object[,hit1:=(Yt==0)*Rtw]
    object[,hit2:=(Yt==1)*Rtw]
    BS <- object[,mean(Rtw)]
    ## FIXME: make sure that sindex cannot be 0
    Int0tdMCsurEffARisk <- MatInt0TcidhatMCksurEff[sindex(object[["time"]],t),]
    ## object[,IC:=object[["hit1"]] + object[["hit2"]] - BS + mean(object[["hit1"]])*Int0tdMCsurEffARisk + colMeans(MatInt0TcidhatMCksurEff*object[["hit2"]])]
    object[,IC:=hit1+hit2-BS + mean(hit1)*Int0tdMCsurEffARisk + colMeans(MatInt0TcidhatMCksurEff*hit2)]
    ## print(list(Int0tdMCsurEffARisk=Int0tdMCsurEffARisk,MatInt0TcidhatMCksurEff=MatInt0TcidhatMCksurEff))
    ## print(object)
    ## [,data.table(hit1,hit2,IC,BS=mean(Rtw))])
    #compute mean and sd of iid representation
    ## sd <- sd(IC)/sqrt(n)
    ## out <- list(epsilon,sd=sd)
    ## sd 
}


getInfluenceCurve.KM <- function(time,status){
    ## compute influence function for reverse Kaplan-Meier
    ## i = 1,..., n are columns
    ## s = 1,..., s_tmax are rows
    N <- length(time)
    times <- unique(time)
    NU <- length(times)
    lagtime <- c(0,times[-NU])
    F <- prodlim(Hist(time,status)~1,data=data.frame(time=time,status=status),reverse=FALSE)
    G <- prodlim(Hist(time,status)~1,data=data.frame(time=time,status=status),reverse=TRUE)
    Stilde.T <- predictSurvIndividual(F,lag=1)*predictSurvIndividual(G,lag=1)
    Stilde.s <- predict(F,times=lagtime)*predict(G,times=lagtime)
    out <- lapply(1:N,function(i){
                      ((1-status[i])*(time[i]<=times))/Stilde.T[i] - cumsum((time[i]>=times)*(G$hazard*G$n.lost)/Stilde.s)
                  })
    do.call("cbind",out)
}

