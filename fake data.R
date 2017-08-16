rm(list=ls(all=TRUE))
set.seed(20)

setwd('U:\\independent studies\\GP approxim')
dat=read.csv('ext covariates.csv',as.is=T)
nobs=nrow(dat)

#get xmat
ind=which(colnames(dat)=='microsc1')
tmp=unique(dat[,-ind]); dim(tmp); max(dat$loc.id)
nloc=nrow(tmp)
ind=which(colnames(tmp)=='loc.id')
xmat=data.matrix(tmp[order(tmp$loc.id),-ind])
p=ncol(xmat)

beta=rnorm(p,mean=0,sd=0.3)
ind=sample(1:p,size=p*0.9)
beta[ind]=0
beta.true=beta

#generate z's
media=xmat[dat$loc.id,]%*%beta
hist(media)
z.true=z=rnorm(nobs,mean=media,sd=1)
dat$microsc1=ifelse(z>0,1,0)

# k=data.frame(z=z.true,mic=dat$microsc1)
# boxplot(z.true~mic,data=k)

write.csv(dat,'fake data.csv',row.names=F)
prev=data.frame(true.prev=pnorm(xmat%*%beta),loc.id=1:nloc)
write.csv(prev,'true prevalences.csv',row.names=F)
