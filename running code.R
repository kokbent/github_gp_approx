rm(list=ls(all=TRUE))
library('mvtnorm')
set.seed(18)

setwd('U:\\independent studies\\GP approxim\\github_gp_approx')
source('gibbs sampler.R')
source('DVfunctions_gibbs.R')

dat=read.csv('fake data.csv',as.is=T)
ngibbs=10000
ind=which(colnames(dat)%in%c('microsc1','loc.id','interc'))
nomes.cov=colnames(dat)[-ind]

#fit model
model.fit=gibbs.approx(dat=dat,nomes.cov=nomes.cov,ngibbs=ngibbs)  
