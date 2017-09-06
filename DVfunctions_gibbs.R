tnorm <- function(n,lo,hi,mu,sig){   #generates truncated normal variates based on cumulative normal distribution
  #normal truncated lo and hi
  
  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))
  
  q1 <- pnorm(lo,mu,sig) #cumulative distribution
  q2 <- pnorm(hi,mu,sig) #cumulative distribution
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  z[z == -Inf]  <- lo[z == -Inf]
  z[z == Inf]   <- hi[z == Inf]
  z
}
#------------------------------------
rmvnorm=function (n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), 
                  method = c("eigen", "svd", "chol")) 
{
  if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
                   check.attributes = FALSE)) {
    stop("sigma must be a symmetric matrix")
  }
  if (length(mean) != nrow(sigma)) {
    stop("mean and sigma have non-conforming size")
  }
  sigma1 <- sigma
  dimnames(sigma1) <- NULL
  if (!isTRUE(all.equal(sigma1, t(sigma1)))) {
    warning("sigma is numerically not symmetric")
  }
  method <- match.arg(method)
  if (method == "eigen") {
    ev <- eigen(sigma, symmetric = TRUE)
    if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
      warning("sigma is numerically not positive definite")
    }
    retval <- ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*% 
      t(ev$vectors)
  }
  else if (method == "svd") {
    sigsvd <- svd(sigma)
    if (!all(sigsvd$d >= -sqrt(.Machine$double.eps) * abs(sigsvd$d[1]))) {
      warning("sigma is numerically not positive definite")
    }
    retval <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
  }
  else if (method == "chol") {
    retval <- chol(sigma, pivot = TRUE)
    o <- order(attr(retval, "pivot"))
    retval <- retval[, o]
  }
  retval <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*% retval
  retval <- sweep(retval, 2, mean, "+")
  colnames(retval) <- names(mean)
  retval
}
#------------------------------------
update.betas=function(cov,lambda,z){
  p=ncol(cov)
  Tinv=diag(x=1/lambda,p)
  Tinv[1,1]=1/10
  prec=t(cov)%*%cov+Tinv
  var1=solve(prec)
  pmedia=t(cov)%*%z
  rmvnorm(1,var1%*%pmedia,var1)
}
#------------------------------------
update.z=function(cov,betas,dat,nobs){
  cond=dat$microsc1>0
  media=cov%*%betas
  res=rep(NA,nobs)
  res[cond]=tnorm(sum(cond),lo=0,hi=Inf,mu=media[cond],sig=1)
  res[!cond]=tnorm(sum(!cond),lo=-Inf,hi=0,mu=media[!cond],sig=1)
  res
}
#------------------------------------
update.lambda=function(cov,betas,a.lamb,b.lamb){
  p=ncol(cov)
  a1=a.lamb+(p-1)/2
  b1=b.lamb+sum(betas[-1]^2)/2
  1/rgamma(1,a1,b1)
}
#------------------------------------
log.marg.likel=function(cov,z,lambda){
  p=ncol(cov)
  
  Tinv=diag(x=1/lambda,p)
  Tinv[1,1]=1/10
  prec=t(cov)%*%cov+Tinv
  var1=solve(prec)
  mu=var1%*%t(cov)%*%z
  
  diag1=c(10,rep(lambda,p-1))
  (1/2)*(-sum(log(diag1))-t(z)%*%z+t(mu)%*%prec%*%mu+determinant(var1)$modulus[1])
}
#------------------------------------
samp.move=function(indin,indout,maxp,cov,z,lambda,xmat.orig){
  indin.old=indin
  p=length(indin.old)
  rand1=runif(1)	
  p0=1
  if (p == 1) {
    indin.new=birth(indin,indout)
    p0=1/3 #death prob 2 -> 1 is (1/3) and birth prob 1 -> 2 is 1. 
  }
  if (p == maxp) {
    if (rand1 < 1/2) {
      indin.new=death(indin)
      p0=2/3 #birth prob T-1 -> T is (1/3) and death prob T -> T-1 is 1/2
    }
    if (rand1 >= 1/2) indin.new=swap(indin,indout)
  }
  if (1 < p & p < maxp) {
    if (rand1 < 1/3) {
      indin.new=birth(indin,indout)
      if (p==maxp-1) p0=3/2 #death prob from T -> T-1 is (1/2) and birth prob from T-1 -> T is (1/3)
    }
    if (1/3 < rand1 & rand1 < 2/3) {
      indin.new=death(indin)
      if (p==2) p0=3 #birth prob from 1 -> 2 is 1 and death prob from 2 -> 1 is 1/3
    }
    if (2/3 < rand1) indin.new=swap(indin,indout)
  }
  pold=log.marg.likel(cov=xmat.orig[,indin.old],z=z,lambda=lambda)
  pnew=log.marg.likel(cov=xmat.orig[,indin.new],z=z,lambda=lambda)+log(p0)
  prob=exp(pnew-pold)
  rand2=runif(1)
  
  seq1=1:maxp
  k=which(!seq1%in%indin.new)
  indout.new=seq1[k]
  if (rand2<prob) return(list(xmat=xmat.orig[,indin.new],indin=indin.new,indout=indout.new))
  return(list(xmat=xmat.orig[,indin.old],indin=indin.old,indout=indout))
}
#------------------------------------------
death=function(indinz){
  k=sample(2:length(indinz),size=1) 
  indinz[-k]
}
#---------------------------------------------------------------------------------------------------
swap=function(indinz,indoutz){
  if (length(indinz)==2) k=indinz[2] #cannot swap intercept
  if (length(indinz)!=2) k=sample(2:length(indinz),size=1)  
  tmp=indinz[-k]
  include=sample(indoutz,size=1)
  sort(c(tmp,include))
}
#---------------------------------------------------------------------------------------------------
birth=function(indinz,indoutz){
  k=sample(indoutz,size=1)
  sort(c(indinz,k))
}
