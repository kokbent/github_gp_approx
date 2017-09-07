library(stringr)

GLMEnvLL <- function (mTrain, mHold, lHold, mod) {
  modTrain <- glm(mod$formula, data = mTrain, family = binomial(link = "probit"))
  n <- as.numeric(table(mHold$cluster))
  m <- as.numeric(table(mHold$cluster, mHold$micro_mala)[,2])
  
  if (any(str_detect(deparse(mod$formula), "age"))) {
    ageRange <- length(unique(mTrain$age))
    seqToRep <- rep(1:nrow(lHold), each = ageRange)
    lHoldExp <- lHold[seqToRep,]
    lHoldExp$age <- rep(c(0:(ageRange-1)), nrow(lHold))
    
    ageProbBase <- 1/(ageRange * 2 - 1)
    ageProb <- ifelse(lHoldExp$age == 0, ageProbBase, 2 * ageProbBase)
    
    prob <- predict(modTrain, lHoldExp, type = "response")
    prob <- prob * ageProb
    prob <- aggregate(prob ~ lHoldExp$cluster, FUN = "sum") %>%
      .[,2]
  }
  else {
    prob <- predict(modTrain, lHold, type = "response")
  }
  
  LL <- sum(dbinom(m, n, prob, log = T))
  return(LL)
}

calcSocProbMN <- function (trainDat, holdDat, wealDecay, eduDecay, eduCat = 0:3) {
  wealMod <- multinom(wFormula,
                      data = trainDat, decay = wealDecay, maxit = 500, trace = F)
  wealHoldProb <- predict(wealMod, holdDat, "probs")
  
  ldatHoldGrid1 <- expand.grid(cluster = holdDat$cluster, wealth = 1:5) %>%
    inner_join(holdDat) %>%
    arrange(cluster, wealth)
  ldatHoldGrid1$wealthProb <- as.vector(t(wealHoldProb))
  
  eduMod <- multinom(eduFormula,
                     data = trainDat, decay = eduDecay, maxit = 500, trace = F)
  
  eduHoldProb <- predict(eduMod, ldatHoldGrid1, "probs")
  
  ldatHoldGrid2 <- expand.grid(cluster = holdDat$cluster, Edu = eduCat) %>%
    inner_join(ldatHoldGrid1) %>%
    arrange(cluster, wealth, Edu)
  ldatHoldGrid2$eduProb <- as.vector(t(eduHoldProb))
  
  return(ldatHoldGrid2)
}

calcFullProb <- function (grid, margProb) {
  grid$fullProb <- with(grid, eduProb * wealthProb * margProb)
  
  fullProb <- aggregate(fullProb ~ cluster, data = grid, FUN = sum) %>%
    .$fullProb
  
  return(fullProb)
}

calcSocProbNN <- function (trainDat, holdDat, wealDecay, eduDecay, eduCat = 0:3) {
  wealMod <- avNNet(wFormula,
                    data = trainDat, repeats = 10, decay = wealDecay, size = 2, 
                    maxit = 500, trace = F, allowParallel = T)
  
  wealHoldProb <- predict(wealMod, holdDat, "prob")
  
  ldatHoldGrid1 <- expand.grid(cluster = holdDat$cluster, wealth = 1:5) %>%
    inner_join(holdDat) %>%
    arrange(cluster, wealth)
  ldatHoldGrid1$wealthProb <- as.vector(t(wealHoldProb))
  
  eduMod <- avNNet(eduFormula,
                   data = trainDat, repeats = 10, decay = eduDecay, size = 2, 
                   maxit = 500, trace = F, allowParallel = T)
  
  eduHoldProb <- predict(eduMod, ldatHoldGrid1, "prob")
  
  ldatHoldGrid2 <- expand.grid(cluster = holdDat$cluster, Edu = eduCat) %>%
    inner_join(ldatHoldGrid1) %>%
    arrange(cluster, wealth, Edu)
  ldatHoldGrid2$eduProb <- as.vector(t(eduHoldProb))
  
  return(ldatHoldGrid2)
}

calcSocProbMixed <- function (trainDat, holdDat, wealDecay, eduDecay, eduCat = 0:3) {
  wealMod <- avNNet(wFormula,
                    data = trainDat, repeats = 10, decay = wealDecay, size = 2, 
                    maxit = 500, trace = F, allowParallel = T)
  
  wealHoldProb <- predict(wealMod, holdDat, "prob")
  
  ldatHoldGrid1 <- expand.grid(cluster = holdDat$cluster, wealth = 1:5) %>%
    inner_join(holdDat) %>%
    arrange(cluster, wealth)
  ldatHoldGrid1$wealthProb <- as.vector(t(wealHoldProb))
  
  eduMod <- polr(eduFormula2, data = trainDat)
  
  eduHoldProb <- predict(eduMod, ldatHoldGrid1, "prob")
  
  ldatHoldGrid2 <- expand.grid(cluster = holdDat$cluster, Edu = eduCat) %>%
    inner_join(ldatHoldGrid1) %>%
    arrange(cluster, wealth, Edu)
  ldatHoldGrid2$eduProb <- as.vector(t(eduHoldProb))
  
  return(ldatHoldGrid2)
}

calcSocProbSimpMixed <- function (trainDat, holdDat, wealDecay, eduDecay, eduCat = 0:3) {
  wealMod <- multinom(wFormula,
                      data = trainDat, decay = wealDecay, maxit = 500, trace = F)
  wealHoldProb <- predict(wealMod, holdDat, "probs")
  
  ldatHoldGrid1 <- expand.grid(cluster = holdDat$cluster, wealth = 1:5) %>%
    inner_join(holdDat) %>%
    arrange(cluster, wealth)
  ldatHoldGrid1$wealthProb <- as.vector(t(wealHoldProb))
  
  eduMod <- polr(eduFormula, data = trainDat)
  
  eduHoldProb <- predict(eduMod, ldatHoldGrid1, "prob")
  
  ldatHoldGrid2 <- expand.grid(cluster = holdDat$cluster, Edu = eduCat) %>%
    inner_join(ldatHoldGrid1) %>%
    arrange(cluster, wealth, Edu)
  ldatHoldGrid2$eduProb <- as.vector(t(eduHoldProb))
  
  return(ldatHoldGrid2)
}

classifyAge <- function (v) {
  age <- v %/% 12
  return(age)
}

calcFullProbAge <- function (grid, margProb) {
  grid$fullProb <- with(grid, eduProb * wealthProb * ageProb * margProb)
  
  fullProb <- aggregate(fullProb ~ cluster, data = grid, FUN = sum) %>%
    .$fullProb
  
  return(fullProb)
}

gridAddAge <- function (grid) {
  gridAge <- expand.grid(cluster = sort(unique(grid$cluster)), age = 0:4) %>%
    inner_join(grid) %>%
    arrange(cluster, wealth, Edu, age)
  
  ind <- gridAge$age == 0
  gridAge$ageProb <- NA
  gridAge$ageProb[ind] <- 1/9
  gridAge$ageProb[!ind] <- 2/9
  
  return(gridAge)
}

gridAddAgeGH <- function (grid) {
  gridAge <- expand.grid(cluster = sort(unique(grid$cluster)), age = 0:5) %>%
    inner_join(grid) %>%
    arrange(cluster, wealth, Edu, age)
  
  gridAge$ageProb <- NA
  ind1 <- gridAge$age == 0
  ind2 <- gridAge$age == 5
  gridAge$ageProb[ind1] <- 1/11
  gridAge$ageProb[ind2] <- 2/11
  gridAge$ageProb[!ind1 & !ind2] <- 2/11
  
  return(gridAge)
}

calcStdLL <- function (actual, modelProb, param) {
  p <- actual[,c(param)]
  actualMat <- as.matrix(table(actual$cluster, p))
  
  dm <- 0
  for (z in 1:nrow(actualMat)) {
    dm <- sum(dm, dmultinom(actualMat[z,], prob = modelProb[z,], log = T))
  }
  dm <- dm/sum(actualMat)
  return(dm)
}

nCV10 <- function (dat, shuff, h) {
  n <- rep(NA, 10)
  for (t in 1:10)  {
    ord <- shuff[((t-1)*h+1):(t*h)] %>%
      sort()
    
    tmp <- dat$cluster %in% ord
    
    n[t] <- sum(tmp)
  }
  
  return(n)
}

gridAddAge2 <- function (grid) {
  gridAge <- expand.grid(cluster = sort(unique(grid$cluster)), age = 0:4) %>%
    inner_join(grid) %>%
    arrange(cluster, wealth, momEdu, age)
  
  ind <- gridAge$age == 0
  gridAge$ageProb <- NA
  gridAge$ageProb[ind] <- 1/9
  gridAge$ageProb[!ind] <- 2/9
  
  return(gridAge)
}
