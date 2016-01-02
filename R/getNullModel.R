getNullModel <- function(formula,data,responseType){
  nullform <- reformulate("1",response=formula[[2]])
  nullfit <- switch(responseType,
                    "binary"={
                      list("Prevalence"=glm(nullform,data,family="binomial"))
                    },
                    "continuous"={
                      stop("FIXME: Need to adapt a function for the empirical distribution function or glm.")
                    },
                    "survival"={
                      list("Kaplan-Meier"=prodlim(nullform,data))
                    },
                    "competing.risks"={
                      list("AalenJohansen"=prodlim(nullform,data))
                    })
  nullfit[[1]]$call$data <- NULL
  nullfit[[1]]$formula <- NULL
  nullfit[[1]]$call$formula=nullform
  nullfit
}
  
