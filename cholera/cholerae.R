library(ape)
library(mlesky)
rm(list=ls())
set.seed(0)
t=read.tree('cholerae.nwk')
t$root.time=1950

res=optim_res_aic(t,ncpu=6,res=10:50,model=2)
print(res)
fit=mlskygrid(t,res=res,tau=NULL,tau_lower = 0.001,tau_upper = 10000,model = 2,ncpu=6)

pdf('cholerae.pdf',7,8)
par(mfrow=c(2,1),mar=c(4,4,1,4))
plot(t,show.tip.label=F)
axisPhylo(1,backward = F)
plot(fit,logy = F)
dev.off()