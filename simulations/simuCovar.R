library(ape)
library(mlesky)
rm(list=ls())

betas=sigmas=rep(NA,100)
for (b in seq(1,length(betas),1)) {
  set.seed(b)
  print(b)
  x=seq(1990,2020,1/12)#monthly
  driver=-((x-2005))^2/100+1
  driver[1:floor(length(driver)/2)]=pmax(driver[1:floor(length(driver)/2)],0)
  sigma=floor((b-1)/10)
  #rho=driver+rnorm(length(driver),0,sigma)
  rho=driver*sigma
  ne=rep(1,length(x))
  for (i in 2:length(x)) ne[i]=ne[i-1]*(1+rho[i-1]*(x[2]-x[1]))
  
  alphaFun=function(newx){approx(x,ne,newx,rule=2)$y}
  
  sampleDates=seq(2000,2020,0.1)
  tree=simCoal(sampleDates,alphaFun,alphaMin = max(0.1,min(ne)))
  
  if (b==11) {
    pdf('simuCovar.pdf',7,10)
    par(mfrow=c(3,1))
    plot(x,driver,type = 'l',xlab = '',ylab='Covariate')
#    plot(x,rho,type = 'l',xlab = '',ylab='Growth rate')
    plot(x,ne,type='l',xlab = '',ylab='Effective population size')
    plot(tree,show.tip.label = F)
    axisPhylo(1,backward = F)
    dev.off()
  }
  
  covar=cbind(x,driver)
  colnames(covar)<-c('time','var')
  covar=as.data.frame(covar)
  sampleTimes=tree$root.time+dist.nodes(tree)[Ntip(tree)+1,]
  fit=mlskygrid(tree,res=30,tau=10,model = 2,ncpu=6,sampleTimes=sampleTimes,formula=diffLogNe~var,data=covar)
  betas[b]=fit$beta[2]
  sigmas[b]=sigma
}

pdf('simuCovar2.pdf',10,10)
boxplot(betas~sigmas,ylab='Inferred association coefficient (beta)',xlab='Multiplicative factor between covariate data and growth rate',outline=F)
dev.off()

#par(mfrow=c(3,1))
#plot(fit$time,fit$ne,type='l')
#plot(fit$time,fit$growthrate,type='l')
#plot(x,driver,type='l')
#print(fit$beta)
