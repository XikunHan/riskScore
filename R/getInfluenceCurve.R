## IMPORTANT : data have to be ordered by time (and for ties by reverse status)
getInfluenceCurve.AUC.survival <- function(t,n,time,status,risk,times,ipcwControls1,ipcwCases,Cases,Controls1,MatInt0TcidhatMCksurEff){
    F01t <- sum(ipcwCases)
    St <- sum(ipcwControls1)
    nbCases <- sum(Cases)
    nbControls1 <- sum(Controls1)
    mcase <- matrix(risk[Cases],nrow=nbCases,ncol=nbControls1)
    mcontrol <- matrix(risk[Controls1],nrow=nbCases,ncol=nbControls1,byrow=TRUE)
    wcase <- matrix(ipcwCases[Cases],nrow=nbCases,ncol=nbControls1)
    wcontrol <- matrix(ipcwControls1[Controls1],nrow=nbCases,ncol=nbControls1,byrow=TRUE)
    Mathtij1 <- (1*(mcase>mcontrol)+.5*(mcase==mcontrol))*wcase*wcontrol*n*n
    ht <- (sum(Mathtij1))/(n*n) 
    vectdit <- Cases*ipcwCases*n
    colSumsMathtij1 <- rep(0,n) # initialise at 0
    colSumsMathtij1[Cases] <- rowSums(Mathtij1) 
    rowSumsMathtij1 <- rep(0,n) # initialize at 0
    rowSumsMathtij1[Controls1] <- colSums(Mathtij1)
    hathtstar <- (sum(Mathtij1))/(n*n)  
    vectTisupt <- n*Controls1/sum(Controls1)
    T1 <- colSums(crossprod(Mathtij1,1+MatInt0TcidhatMCksurEff[Cases,]))/n
    T3 <- colSums(hathtstar*(vectTisupt + (vectdit*(1+MatInt0TcidhatMCksurEff)-F01t)/F01t))
    Term.ijak <- (T1-T3)/(F01t*St)
    Term.ikaj <- (rowSumsMathtij1 - n*hathtstar)/(F01t*St)
    Term.jkai <-  (colSumsMathtij1 - n*hathtstar*(vectTisupt+(1/F01t)*(vectdit-F01t)))/(F01t*St)
    ## the influence function according to Blanche et al. 2013, DOI: 10.1002/sim.5958, Statistics in Medicine, Appendix A
    (Term.ijak + Term.ikaj + Term.jkai)/(n)
}

getInfluenceCurve.Brier <- function(t,time,Yt,ipcwResiduals,MatInt0TcidhatMCksurEff){
    browser()
    hit1=(Yt==0)*ipcwResiduals
    hit2=(Yt==1)*ipcwResiduals
    Brier <- mean(ipcwResiduals)
    ## FIXME: make sure that sindex cannot be 0
    Int0tdMCsurEffARisk <- MatInt0TcidhatMCksurEff[prodlim::sindex(time,t),]
    IC.Brier=hit1+hit2-Brier + mean(hit1)*Int0tdMCsurEffARisk + colMeans(MatInt0TcidhatMCksurEff*hit2)
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

