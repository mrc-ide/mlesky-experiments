library(ape)
library(mlesky)
rm(list=ls())
set.seed(2)
alphaFun=function(x){if (x<2005 || x>2010) 10 else 1}
sampleDates=seq(2000,2020,0.1)
tree=simCoal(sampleDates,alphaFun,alphaMin = 0.1)
sts=tree$root.time+dist.nodes(tree)[Ntip(tree)+1,1:Ntip(tree)]
names(sts)=tree$tip.label

res=suggest_res(tree)#optim_res_aic(tree,ncpu=6,res=1:50,model=2)
print(res)
fit=mlskygrid(tree,sampleTimes = sts,res=res,tau=NULL,tau_lower = 0.001,tau_upper = 10000,model = 2,ncpu=6)
fit=parboot(fit,nrep=100,ncpu=6)

pdf('simuBottle.pdf',7,10)
par(mfrow=c(2,1),mar=c(4,4,1,4))
if (!is.null(tree$root.time)) from=tree$root.time else from=-max(dist.nodes(tree)[Ntip(tree)+1,])
to=from+max(dist.nodes(tree)[Ntip(tree)+1,])
xs=seq(from=from,to=to,length.out=100)
ys=xs
for (i in 1:length(ys)) ys[i]=alphaFun(ys[i])
plot(tree,show.tip.label = F)
axisPhylo(1,backward = F)
plot(fit)
lines(xs,ys,col='red')
dev.off()
system('open simuBottle.pdf')

#distance=function(alphaFun,fit) {
#  mean((fit$ne-alphaFun(fit$time+max(sampleDates)))^2)
#}
#print(distance(alphaFun,fit))
