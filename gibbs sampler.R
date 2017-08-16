gibbs.approx=function(dat=dat,ngibbs=ngibbs){
  
nobs=nrow(dat)
ind=which(colnames(dat)%in%c('microsc1','loc.id'))
xmat.orig=cov=data.matrix(dat[,-ind])
maxp=ncol(cov)

#initial values
betas=c(0.2,rep(0,maxp-1))
cond=dat$microsc1>0
z=ifelse(cond,runif(nobs),runif(nobs,min=-1,max=0))
indin=1:10
indout=11:maxp
betas=betas[indin]
cov=cov[,indin]
lambda=0.1

#prior
a.lamb=0.999; b.lamb=0.00999
# a.lamb=2.349; b.lamb=0.2349
# a.lamb=0.1; b.lamb=0.1

ngibbs=10000
vec.betas=matrix(NA,ngibbs,maxp)
vec.outros=matrix(NA,ngibbs,1)

for (i in 1:ngibbs){
  print(c(i,indin))
  if (!1%in%indin) break;
  tmp=samp.move(indin=indin,indout=indout,maxp=maxp,cov=cov,z=z,
                lambda=lambda,xmat.orig=xmat.orig)
  cov=tmp$xmat
  indin=tmp$indin
  indout=tmp$indout
  
  betas=t(update.betas(cov=cov,lambda=lambda,z=z))
  
  lambda=update.lambda(cov=cov,betas=betas,a.lamb=a.lamb,b.lamb=b.lamb)
  z=update.z(cov=cov,betas=betas,dat=dat,nobs=nobs)
  
  tmp=rep(0,maxp)
  tmp[indin]=betas
  vec.betas[i,]=tmp
  vec.outros[i]=lambda
}

list(betas=vec.betas,lambda=vec.outros)
# seq1=500:(i-1)
# plot(vec.betas[seq1,1],type='l')
# plot(vec.outros[seq1,1],type='l')
# plot(apply(vec.betas[seq1,]==0,2,mean),type='h')

}