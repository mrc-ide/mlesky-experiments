library(ape)
library(mlesky)
rm(list=ls())
set.seed(1)
alphaFun=function(x){20}
sampleDates=seq(2000,2020,0.1)
t=simCoal(sampleDates,alphaFun,alphaMin = 0.1)
pdf('simuConst.pdf',7,7)
plotBoth(t,alphaFun)
dev.off()

res=optim_res_aic(t,ncpu=6,res=1:20,model=2)
print(res)

pdf('simuConstResult.pdf',7,10)
par(mfrow=c(3,3))
for (res in c(5,20,50)) for (tau in c(1,10,20)) {
fit=mlskygrid(t,res=res,tau=tau,model = 2,ncpu=6)
plot(fit,logy = F,ylim=c(0,100))
}
dev.off()
