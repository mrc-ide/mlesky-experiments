library(ape)
library(mlesky)
rm(list=ls())
set.seed(1)
alphaFun=function(x){20}
sampleDates=seq(2000,2020,0.1)
res=c()
for (i in 1:100) {
t=simCoal(sampleDates,alphaFun,alphaMin = 0.1)
#res=c(res,suggest_res(t))
res=c(res,optim_res_aic(t,ncpu=6,res=1:20,model=2))
}
length(which(res==1))
