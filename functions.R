create_K <- function (l, n) {
  N <- sum(n)
  K <- matrix(0, nrow = N, ncol = l)
  for (i in 1:l) {
    ini <- sum(n[0:(i-1)]) + 1
    fin <- sum(n[0:i])
    K[ini:fin,i] <- 1
  }
  return(K)
}

create_Sigma <- function (r, D) {
  S_r <- -D/r
  S_r <- exp(S_r)
  return(S_r)
}

tnorm <- function(n,lo,hi,mu,sig){   #generates truncated normal variates based on cumulative normal distribution
  #normal truncated lo and hi
  
  if(length(lo) == 1 & length(mu) > 1) lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1) hi <- rep(hi,length(mu))
  
  q1 <- pnorm(lo,mu,sig) #cumulative distribution
  q2 <- pnorm(hi,mu,sig) #cumulative distribution
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  z[z == -Inf]  <- lo[z == -Inf]
  z[z == Inf]   <- hi[z == Inf]
  z
}

update_beta <- function (param) {
  D10 <- diag(10, nrow = p_x)
  Psi_b <- solve(t(X) %*% X +  solve(D10))
  mu_b <- Psi_b %*% (t(X) %*% (param$z - K %*% param$a))
  return(as.vector(rmvnorm(1, mu_b, Psi_b)))
}

update_gamma <- function (p, s) {
  Psi_g <- solve(t(s$E) %*% ((1/p$sigmasq) * p$Sigma_g_inv) %*% s$E + s$D10_inv)
  mu_g <-  Psi_g %*% (t(s$E) %*% ((1/p$sigmasq) * p$Sigma_g_inv) %*% p$a)
  return(as.vector(rmvnorm(1, mu_g, Psi_g)))
}

update_a <- function (p, s) {
  Psi_a <- solve(s$KK_mat + (1/p$sigmasq) * p$Sigma_g_inv)
  mu_a <- Psi_a %*% (t(s$K) %*% (p$z) + 
                       (1/p$sigmasq) * p$Sigma_g_inv %*% s$E %*% p$gamma)
  return(as.vector(rmvnorm(1, mu_a, Psi_a)))
}

update_z <- function(p, s) {
  mean_z <- s$K %*% p$a
  z.lo <- ifelse(s$m>0, 0, -Inf)
  z.hi <- ifelse(s$m>0, Inf, 0)
  z <- tnorm(s$N, z.lo, z.hi, mean_z, 1)
  
  if (sum(z == Inf | z == -Inf) > 0) print("panic!")
  z[z==Inf] <- 10
  z[z==-Inf] <- -10
  return(as.vector(z))
}

rho_likelihood <- function (r) {
  Sigma <- create_Sigma(r)
  mu <- params_g$a - E %*% params_g$gamma
  l <- log(det(Sigma)^(-0.5)) + (-0.5 * t(mu) %*% solve(Sigma) %*% mu)
  return(l)
}

proposal_likelihood <- function (x, mu, sig) {
  # tnorm(0, Inf)
  lo <- pnorm(0, mu, sig)
  hi <- pnorm(Inf, mu, sig)
  L <- dunif(pnorm(x, mu, sig), lo, hi)
  return(L)
}

update_rho <- function(param) {
  # Metropolis algorithm (No Hasting yet)
  rho_new <- param$rho
  rho_proposed <- tnorm(1, 0, Inf, param$rho, 2)
  lhr <- rho_likelihood(rho_proposed) - rho_likelihood(param$rho) +
    log(proposal_likelihood(param$rho, rho_proposed, 2)) -
    log(proposal_likelihood(rho_proposed, param$rho, 2))
  
  if (log(runif(1)) < lhr) {
    rho_new <- rho_proposed
    rho_accept <<- rho_accept+1
  }
  
  return(rho_new)
}

update_rho_mn <- function(p, s) {
  err <- p$a - s$E %*% p$gamma
  
  like <- rep(NA, length(s$rho_cand))
  for (i in 1:length(s$rho_cand)) {
    like[i] <- t(err) %*% s$Sigma_inv_cand[[i]] %*% err
  }
  
  like <- -0.5 * s$lndetS_cand - ((s$l-1)/2)*log(like)
  weight <- exp(like - max(like))
  weight <- weight/sum(weight)
  
  tmp <- rmultinom(1, 1, weight)
  k <- which(tmp == 1)
  
  return(k)
}

update_sigmasq <- function(p, s) {
  # alpha <- (l + 0.02)/2
  # beta <- 0.5 * t(p$a - E %*% p$gamma) %*% Sigma_g_inv %*% (p$a - E %*% p$gamma) + 0.02
  alpha <- (s$l-1)/2
  err <- p$a - s$E %*% p$gamma
  beta <- 0.5 * t(err) %*% p$Sigma_g_inv %*% err
  return(1/rgamma(1, alpha, rate=beta))
}

library(tcltk)

fragmentDat <- function (clusters, wdat, ldat, latlong, mdat) {
  ord <- sort(clusters)
  datList <- list(wdat = list(), ldat = list(), latlong = list(), mdat = list())
  
  tmp <- wdat$cluster %in% ord
  datList$wdat$hold <- wdat[tmp,]
  datList$wdat$train <- wdat[!tmp,]
  
  tmp <- ldat2$cluster %in% ord
  datList$ldat$hold <- ldat[tmp,]
  datList$ldat$train <- ldat[!tmp,]
  
  tmp <- latlong$cluster %in% ord
  datList$latlong$hold <- latlong[tmp,]
  datList$latlong$train <- latlong[!tmp,]
  
  tmp <- mdat$cluster %in% ord
  datList$mdat$hold <- mdat[tmp,]
  datList$mdat$train <- mdat[!tmp,]
  
  return(datList)
}

gibbsEnv <- function(datList, call, iter, betaPrior = 0.1) {
  # Setting up for Gibbs sampler
  l <- nrow(datList$ldat$train)
  n <- as.vector(table(datList$mdat$train$cluster))
  N <- sum(n)
  
  E <- model.matrix(call, datList$ldat$train)
  m <- datList$mdat$train$micro_mala
  p_e <- ncol(E)
  
  K <- create_K(l, n)
  KK_mat <- t(K) %*% K
  D <- datList$latlong$train[, c("LONGNUM", "LATNUM")]
  D_mat <- as.matrix(dist(D))
  d <- c(10, rep(betaPrior, p_e - 1))
  D10_inv <- diag(1 / d)
  
  # Discretizing rho
  cor1 <- seq(from = 0.01, to = 0.90, by = 0.01)
  # rho_cand <- seq(from = 0.1, to = 10, by = 0.1)
  rho_cand <- -quantile(D_mat, 0.1) / log(cor1)
  Sigma_cand <- list()
  Sigma_inv_cand <- list()
  lndetS_cand <- rep(NA, length(rho_cand))
  for (i in 1:length(rho_cand)) {
    Sigma_cand[[i]] <- create_Sigma(rho_cand[i], D_mat)
    lndetS_cand[i] <- determinant(Sigma_cand[[i]])$modulus
    Sigma_inv_cand[[i]] <- solve(Sigma_cand[[i]])
  }
  
  # Gibbs samplers parameter and priors
  rho_g <- sample(rho_cand, 1)
  k <- which(rho_cand == rho_g)
  gamma_g <- as.vector(rmvnorm(1, rep(0, p_e), diag(10, nrow = p_e)))
  sigmasq_g <- (runif(1, 0, 3)) ^ 2
  Sigma_g <- Sigma_cand[[k]]
  Sigma_g_inv <- Sigma_inv_cand[[k]]
  a_g <- as.vector(rmvnorm(1, E %*% gamma_g, sigmasq_g * Sigma_g))
  z_g <- ifelse(m == 0, -1, 1)
  
  # Parameter list
  params_g <- list(gamma = gamma_g,
                   z = z_g,
                   rho = rho_g,
                   sigmasq = sigmasq_g,
                   a = a_g,
                   Sigma = Sigma_g, 
                   Sigma_inv = Sigma_g_inv
  )
  
  # Static list (for passing around)
  staticList <- list(l = l, n = n, N = N, E = E, m = m, K = K, KK_mat = KK_mat,
                     D = D, D_mat = D_mat, D10_inv = D10_inv, rho_cand = rho_cand,
                     Sigma_cand = Sigma_cand, Sigma_inv_cand = Sigma_inv_cand,
                     lndetS_cand = lndetS_cand)
  
  # Storage
  as_g <- matrix(NA, nrow = iter, ncol = l)
  gammas_g <- matrix(NA, nrow = iter, ncol = p_e)
  # zs_g <- matrix(NA, nrow = iter, ncol = N)
  sigmasqs_g <- rep(NA, iter)
  rhos_g <- rep(NA, iter)
  
  for (i in 1:iter) {
    title <- paste("Gibbs Iterations ", t, sep = "")
    if (!exists("pb")) pb <- tkProgressBar(title, min = 1, max = iter)
    setTkProgressBar(pb, i)
    
    ind <- update_rho_mn(params_g, staticList)
    params_g$rho <- rho_cand[ind]
    params_g$Sigma_g <- Sigma_cand[[ind]]
    params_g$Sigma_g_inv <- Sigma_inv_cand[[ind]]
    
    params_g$sigmasq <- update_sigmasq(params_g, staticList)
    params_g$gamma <- update_gamma(params_g, staticList)
    params_g$a <- update_a(params_g, staticList)
    params_g$z <- update_z(params_g, staticList)
    
    gammas_g[i,] <- params_g$gamma
    as_g[i,] <- params_g$a
    rhos_g[i] <- params_g$rho
    # zs_g[i,] <- params_g$z
    sigmasqs_g[i] <- params_g$sigmasq
  }
  close(pb)
  rm(pb)
  
  return(list(a = as_g, gamma =  gammas_g, rho = rhos_g, sigmasq = sigmasqs_g, call = call,
              dataTrain = E))
}

gibbEnvPred <- function (datList, gibbPar) {
  iter <- length(gibbPar$rho)
  burn <- iter / 2
  DTrain <- datList$latlong$train[, c("LONGNUM", "LATNUM")]
  DHold <- datList$latlong$hold[, c("LONGNUM", "LATNUM")]
  D_all <- rbind(DHold, DTrain)
  DAllMat <- as.matrix(dist(D_all))
  l <- nrow(DTrain)
  h <- nrow(DHold)
  lAll <- l + h
  
  E <- gibbPar$dataTrain
  E_exclude <- model.matrix(gibbPar$call, datList$ldat$hold)
  as_pred_g <- matrix(NA, nrow = (iter - burn), ncol = h)
  ms_pred_prob <- matrix(NA, nrow = (iter - burn), ncol = h)
  
  for (i in (burn + 1):iter) {
    title <- paste("Prediction Iterations ", t, sep = "")
    if (!exists("pb1")) pb1 <- tkProgressBar(title, min = 1, max = iter - burn)
    setTkProgressBar(pb1, (i - burn))
    
    Sigma_all <- gibbPar$sigmasq[i] * create_Sigma(gibbPar$rho[i], DAllMat)
    Sigma_12 <- Sigma_all[1:h, (h + 1):lAll]
    Sigma_22 <- Sigma_all[(h + 1):lAll, (h + 1):lAll]
    Sigma_11 <- Sigma_all[1:h, 1:h]
    Sigma1222 <- Sigma_12 %*% solve(Sigma_22)
    mu_pred <- E_exclude %*% gibbPar$gamma[i,] +
      Sigma1222 %*% (as.matrix(gibbPar$a[i,]) - E %*% gibbPar$gamma[i,])
    Sigma_pred <- Sigma_11 - Sigma1222 %*% t(Sigma_12)
    Sigma_pred <- as.matrix(forceSymmetric(Sigma_pred))
    a_pred <- rmvnorm(1, mu_pred, Sigma_pred)
    m_pred <- pnorm(a_pred)
    
    as_pred_g[(i - burn),] <- a_pred
    ms_pred_prob[(i - burn),] <- m_pred
  }
  close(pb1)
  rm(pb1)
  
  return(ms_pred_prob)
}
