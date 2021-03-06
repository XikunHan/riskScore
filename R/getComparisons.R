### getComparisons.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jan  3 2016 (13:30) 
## Version: 
## last-updated: Jan 27 2016 (16:43) 
##           By: Thomas Alexander Gerds
##     Update #: 31
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
getComparisons <- function(dt,NF,N,alpha,dolist=NF:1){
    x=model=IC=NULL
    ## FIXME: when dolist is 0:1 and models are 0:2 this does not work 
    if (length(dolist)>0){
        data.table::rbindlist(lapply(dolist,function(g){
                                         theta <- dt[,list(x=x[1]),by=model]
                                         delta <- theta[model>g][["x"]]-theta[model==g][["x"]]
                                         se.delta <- dt[model>g,list(se=sd(dt[model==g][["IC"]]-IC)/sqrt(N)),by=model][["se"]]
                                         lower <- delta - qnorm(1-alpha/2) * se.delta
                                         upper <- delta + qnorm(1-alpha/2) * se.delta
                                         p <-2*pnorm(abs(delta)/se.delta,lower.tail=FALSE)
                                         data.table(model=theta[model>g][["model"]],reference=g,delta=delta,se.delta=se.delta,lower=lower,upper=upper,p=p)
                                     }))
    }else {NULL}
}
#----------------------------------------------------------------------
### getComparisons.R ends here
