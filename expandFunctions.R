library(stringr)

extractX <- function (str) {
  str_split(str, 'x', simplify = T)
}

addQuadTerm <- function (dat, covLab) {
  colNamesNew <- paste(covLab, "2", sep = "")
  oriDat <- dat[,covLab]
  quadTerms <- oriDat^2
  colnames(quadTerms) <- colNamesNew
  newdat <- cbind(dat, quadTerms)
  return(newdat)
}

addInterTerm <- function(dat, covLab) {
  p <- length(covLab)
  oriDat <- dat[,covLab]
  newdat <- dat
  
  for (i in 1:(p-1)) {
    inter <- oriDat[,i] * oriDat[,-1:-i] %>% as.data.frame()
    colNamesNew <- paste(covLab[i], "x",covLab[-1:-i], sep = "")
    colnames(inter) <- colNamesNew
    newdat <- cbind(newdat, inter)
  }
  return(newdat)
}

addCubiTerm <- function (dat, covLab) {
  colNamesNew <- paste(covLab, "3", sep = "")
  oriDat <- dat[,covLab]
  cubiTerms <- oriDat^3
  colnames(cubiTerms) <- colNamesNew
  newdat <- cbind(dat, cubiTerms)
  return(newdat)
}

addCubiInterTerm <- function(dat, covLab) {
  p <- length(covLab)
  quadName <- paste(covLab, "2", sep = "")
  tmp <- extractX(colnames(dat))
  ind <- tmp[,2] != ""
  inter2Name <- colnames(dat)[ind]
  
  singleTerm <- dat[,covLab]
  quadTerm <- dat[,quadName]
  inter2Term <- dat[,inter2Name]
  
  newdat <- dat
  for (i in 1:p) {
    inter <- singleTerm[,i] * quadTerm[,-i] %>% as.data.frame()
    colNamesNew <- paste(covLab[i], "x",quadName[-i], sep = "")
    colnames(inter) <- colNamesNew
    newdat <- cbind(newdat, inter)
  }
  
  splitInter2Names <- extractX(inter2Name)
  for (i in 1:p) {
    tmp <- splitInter2Names %in% covLab[1:i] %>%
      matrix(ncol=2)
    tmp <- tmp[,1] | tmp[,2]
    if (any(!tmp)) {
      inter <- singleTerm[,i] * inter2Term[,!tmp] %>% as.data.frame()
      colNamesNew <- paste(covLab[i], "x", inter2Name[!tmp], sep = "")
      colnames(inter) <- colNamesNew
      newdat <- cbind(newdat, inter)
    }
  }
  
  return(newdat)
}

removeHiCor <- function (dat, threshold) {
  R <- cor(dat)
  RBool <- abs(R) > threshold & R < 1
  
  hiCorDf <- data.frame(X1 = character(), X2 = character(), stringsAsFactors = F)
  for (i in 1:ncol(RBool)) {
    ind <- RBool[,i]
    if (sum(ind) > 0) {
      X2 <- rownames(RBool)[ind]
      X1 <- rep(colnames(RBool)[i], sum(ind))
      tmpdf <- data.frame(X1 = X1, X2 = X2, stringsAsFactors = F)
      hiCorDf <- rbind(hiCorDf, tmpdf)
    }
  }
  
  i <- 1
  repeat {
    if (i > nrow(hiCorDf)) {
      break
    }
    
    removalCov <- hiCorDf[i,2]
    dat <- dat[,colnames(dat) != removalCov]
    hiCorDf <- hiCorDf[(i+1):nrow(hiCorDf),] %>%
      filter(X1 != removalCov | X2 != removalCov)
    # hiCorDf <- rbind(hiCorDf[1:i,], tmp)
    i <- i+1
  }
  
  return(dat)
}
