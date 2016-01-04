### test-Score.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jan  4 2016 (14:30) 
## Version: 
## last-updated: Jan  4 2016 (15:54) 
##           By: Thomas Alexander Gerds
##     Update #: 4
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
context("riskScore")
test_that("survival outcome: Brier Score",
          {   library(survival)
              library(riskScore)
              library(pec)
              set.seed(112)
              d <- sampleData(112,outcome="survival")
              f1 <- coxph(Surv(time,event)~X1+X5+X8,data=d)
              f2 <- coxph(Surv(time,event)~X2+X6+X9+X10,data=d)
              p1 <- pec(list(f1,f2),formula=Surv(time,event)~1,data=d,times=c(3,5,10),exact=FALSE,start=NULL)
              s1 <- Score(list(f1,f2),formula=Surv(time,event)~1,data=d,times=c(3,5,10),test=FALSE,metrics="brier")
              expect_equal(p1$AppErr$coxph,s1$Brier[model=="coxph",Brier])
              expect_equal(p1$AppErr$coxph.1,s1$Brier[model=="coxph.1",Brier])
              expect_equal(p1$AppErr$Reference,s1$Brier[model=="Kaplan-Meier",Brier])
          })
#----------------------------------------------------------------------
### test-Score.R ends here
