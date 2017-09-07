stopCluster(cl)
rm(list = ls(all = TRUE))

library(mvtnorm)
library(doSNOW)
library(doParallel)
library(iterators)
library(tcltk)
library(Matrix)
library(dplyr)
library(nnet)
cl <- makeCluster(5)
registerDoSNOW(cl)

set.seed(1234)

source("functions.R")
source("cvInitTG.R")
source('gibbs sampler.R')
source('DVfunctions_gibbs.R')

iter <- 20000
fill <- 10 - length(shuff) %% 10
shuffMat <- matrix(c(shuff, rep(NA, fill)), nrow = 10)

predGPA <- function (dat, model.fit) {
  n <- nrow(dat)
  iter <- nrow(model.fit[[1]])
  burn <- iter/2
  seqPred <- (burn+1):iter
  XMat <- as.matrix(dat[,-1])
  XMat <- cbind(1, XMat)
  
  mPredProb <- matrix(NA, nrow = burn, ncol = n)
  for (i in seqPred) {
    muPred <- XMat %*% model.fit$betas[i,]
    aPred <- rnorm(n, muPred, sqrt(model.fit$lambda[i]))
    mPred <- pnorm(aPred)
    mPredProb[(i-burn),] <- mPred
  }
  
  return(mPredProb)
}

ldatExp <- read.csv("Data/TG/TGLocDatExpCubi.csv")
Likelihoods <- foreach(t = 1:10, .packages = c("dplyr", "tcltk"),
                       .errorhandling = "remove") %dopar% 
                       {
                         clust <- shuffMat[t, ] %>% na.exclude()
                         datList <- fragmentDat(clust, wdat, ldat2, latlong, mdat)
                         
                         trainDat <- datList$mdat$train %>% select(cluster, micro_mala) %>%
                           left_join(ldatExp)
                         
                         ngibbs = 10000
                         dat = trainDat
                         colnames(dat)[1] <- 'loc.id'
                         colnames(dat)[2] <- 'microsc1'
                         ind = which(colnames(dat) %in% c('loc.id', 'microsc1'))
                         nomes.cov = colnames(dat)[-ind]
                         
                         model.fit = gibbs.approx(dat = dat,
                                                  nomes.cov = nomes.cov,
                                                  ngibbs = ngibbs)
                         
                         holdDat <- datList$ldat$hold %>% select(cluster) %>%
                           left_join(ldatExp)
                         
                         predProb <- predGPA(holdDat, model.fit)
                         bayesProb <- apply(predProb, 2, mean)
                         
                         nMalaHold <- as.numeric(table(datList$mdat$hold$cluster))
                         mMalaHold <-
                           as.numeric(table(datList$mdat$hold$cluster, datList$mdat$hold$micro_mala)[, 2])
                         
                         LL <- sum(dbinom(mMalaHold, nMalaHold, bayesProb, log = T))
                         LL
                       }
