rm(list=ls(all=TRUE))
library('mvtnorm')
set.seed(18)

setwd('U:\\independent studies\\GP approxim')
source('gibbs sampler.R')
source('DVfunctions_gibbs.R')

dat=read.csv('fake data.csv',as.is=T)
ngibbs=10000

#fit model
model.fit=gibbs.approx(dat=dat,ngibbs=ngibbs)  
