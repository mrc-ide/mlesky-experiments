library(ape)
library(mlesky)
rm(list=ls())

distance=function(alphaFun,fit) {
  mean((fit$ne-alphaFun(fit$time+max(sampleDates)))^2)
}

alphaFun=function(x){sin(x)*10+12}
sampleDates=seq(2000,2020,0.1)
#sampleDates=seq(2000,2020,length.out = 500)

distances=matrix(NA,10,3)
fit=as.list(rep(NA,3))
for (dup in 1:10) {
  set.seed(dup)
  t=simCoal(sampleDates,alphaFun,alphaMin = 0.1)
  for (m in c(1,2,3)) {
    fit[[m]]=mlskygrid(t,res=20,tau=NULL,tau_lower = 0.001,tau_upper = 10000,model = m,ncpu=6)
    distances[dup,m]=distance(alphaFun,fit[[m]])
  }
}

print(c(median(distances[v,2]),median(distances[v,3]),median(distances[v,1])))


