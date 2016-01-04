### getComparisons.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jan  3 2016 (13:30) 
## Version: 
## last-updated: Jan  4 2016 (14:27) 
##           By: Thomas Alexander Gerds
##     Update #: 18
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
getComparisons <- function(dt,NF,N,alpha,dolist=NF:1){
    ## FIXME: when dolist is 0:1 and models are 0:2 this does not work 
    if (length(dolist)>0){
        rbindlist(lapply(dolist,function(g){
                             theta <- dt[,list(x=x[1]),by=model]
                             delta <- theta[model==g][["x"]]-theta[model<g][["x"]]
                             se.delta <- dt[model<g,list(se=sd(dt[model==g][["IC"]]-IC)/sqrt(N)),by=model][["se"]]
                             lower <- delta - qnorm(1-alpha/2) * se.delta
                             upper <- delta + qnorm(1-alpha/2) * se.delta
                             p <-2*pnorm(abs(delta)/se.delta,lower.tail=FALSE)
                             data.table(model1=g,model2=theta[model<g][["model"]],delta=delta,lower=lower,upper=upper,p=p)
                         }))
    }else {NULL}
}


#----------------------------------------------------------------------
### getComparisons.R ends here
