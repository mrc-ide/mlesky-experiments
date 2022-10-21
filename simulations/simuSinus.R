library(ape)
library(mlesky)
rm(list=ls())
set.seed(3)
alphaFun=function(x){sin(x)*10+12}
sampleDates=seq(2000,2020,0.1)
t=simCoal(sampleDates,alphaFun,alphaMin = 0.1)
sts=t$root.time+dist.nodes(t)[Ntip(t)+1,1:Ntip(t)]
names(sts)=t$tip.label
pdf('simu.pdf',7,7)
plotBoth(t,alphaFun)
dev.off()

fit=as.list(rep(NA,3))
for (m in c(1,2,3)) {
  fit[[m]]=mlskygrid(t,sampleTimes = sts,res=20,tau=NULL,tau_lower = 0.001,tau_upper = 10000,model = m,ncpu=6)
  fit[[m]]=parboot(fit[[m]],nrep=100,ncpu=6)
}

pdf('simuResult.pdf',7,7)
par(mfrow=c(3,1))
for (m in c(2,3,1)) {
  plot(fit[[m]],ylim=c(1,100),logy = T)
}
dev.off()
system('open simuResult.pdf')

distance=function(alphaFun,fit) {
  mean((fit$ne-alphaFun(fit$time+max(sampleDates)))^2)
}

distances=c()
for (m in c(2,3,1)) {
  distances=c(distances,distance(alphaFun,fit[[m]]))
}
print(distances)

#fit=mlskygrid(t,res=20,tau=NULL,tau_lower = 0.001,tau_upper = 10000,model = m,ncpu=6)
#fit=parboot(fit,nrep=10,ncpu=6)

