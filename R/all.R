################################
BIC_f_zcut2 <- function(y, rrr, Zcut) {
  n <- length(y)
  nparams <- nrow(rrr) * 3
  w <- rrr[, 1]
  u <- rrr[, 2]
  sig <- rrr[, 3]
  y0 <- y[which(y >= Zcut)]
  cc <- c()
  for (i in 1:nrow(rrr)) {
    c <- dnorm(y0, u[i], sig[i]) * w[i]
    cc <- rbind(cc, c)
  }
  d <- apply(cc, 2, sum)
  e <- sum(log(d))
  f <- e * 2 - nparams * log(n)
  return(f)
}

#' @importFrom LTMGSCA SeparateKRpkmNew
GetBestK <- function(x, n, q, err = 1e-10) {
  best.bic <- -Inf
  best.k <- 0
  best.result <- c(0, 0, 0)
  for (k in 1:7) {
    rrr <- SeparateKRpkmNew(x = x, n = n, q = q, k = k, err = err)
    bic <- BIC_f_zcut2(y = x, rrr, q)
    if (is.nan(bic)) {
      bic <- -Inf
    }
    if (bic >= best.bic) {
      best.bic <- bic
      best.k <- k
      best.result <- rrr
    } else {
      return(list(k = best.k, bic = best.bic, result = best.result))
    }
  }
  return(list(k = 0, bic = 0, result = c(0, 0, 0)))
}

SeparateKRpkmNew2 <- function(x, n = 1000, q = 0, k = 1, err = 1e-10) {
  if (k == 0)
    return(SeparateKRpkmNewSingle(x, n, q, err = 1e-10))
  q <- max(q, min(x))
  c <- sum(x < q)
  x <- x[which(x >= q)]
  if (length(x) <= k) {
    warning(sprintf("The length of x is %i. Sorry, too little conditions\n", length(x)))
    return(cbind(0, 0, 0))
  }
  mean <- c()
  for (i in 1:k) {
    mean <- c(mean, sort(x)[floor(i * length(x)/(k + 1))])
  }
  if (k > 1) {
    mean[1] <- min(x) - 1
    mean[length(mean)] <- max(x) + 1
  }
  p <- rep(1/k, k)
  sd <- rep(sqrt(var(x)), k)
  pdf.x.portion <- matrix(0, length(x), k)
  for (i in 1:n) {
    p0 <- p
    mean0 <- mean
    sd0 <- sd
    pdf.x.all <- t(p0 * vapply(x, function(x) dnorm(x, mean0, sd0), rep(0, k)))
    pdf.x.portion <- pdf.x.all/rowSums(t(pdf.x.all))
    cdf.q <- pnorm(q, mean0, sd0)
    cdf.q.all <- p0 * cdf.q
    cdf.q.portion <- cdf.q.all/sum(cdf.q.all)
    cdf.q.portion.c <- cdf.q.portion * c
    denom <- colSums(t(pdf.x.portion)) + cdf.q.portion.c
    p <- denom/(nrow(t(pdf.x.portion)) + c)
    im <- dnorm(q, mean0, sd0)/cdf.q * sd0
    im[is.na(im)] <- 0
    mean <- colSums(crossprod(x, t(pdf.x.portion)) + (mean0 - sd0 * im) * cdf.q.portion.c)/denom
    sd <- sqrt((colSums((x - matrix(mean0, ncol = length(mean0), nrow = length(x), byrow = TRUE))^2 * t(pdf.x.portion)) + sd0^2 * (1 - (q - mean0)/sd0 * im) *
                  cdf.q.portion.c)/denom)
    if (!is.na(match(NaN, sd))) {
      break
    }
    if ((mean(abs(p - p0)) <= err) && (mean(abs(mean - mean0)) <= err) && (mean(abs(sd - sd0)) <= err)) {
      break
    }
  }
  return(cbind(p, mean, sd))
}

fit_LTMG <- function(x, q, n = 1000, err = 1e-10) {
  best <- GetBestK(x = x, n = n, q = q, err = err)
  return(list(result = best$result, peak.count = best$k))
}

cdf_LTMG <- function(a, ccc) {
  pp <- rep(0, length(a))
  if (sum(abs(ccc)) > 0) {
    for (i in 1:nrow(ccc)) {
      pp <- pp + ccc[i, 1] * pnorm(a, ccc[i, 2], ccc[i, 3])
    }
  }
  return(pp)
}


fit_LTG <- function(a, Zcut0) {
  ccc <- SeparateKRpkmNew2(x = a, q = Zcut0, n = 1000, k = 1, err = 1e-10)
  return(ccc)
}

cdf_LTG <- function(a, ccc) {
  return(pnorm(a, ccc[2], ccc[3]))
}

fit_G <- function(a) {
  mean0 <- mean(a)
  sd0 <- sd(a)
  ccc <- list(mean0, sd0)
  names(ccc) <- c("mean", "sd")
  return(ccc)
}

cdf_G <- function(a, ccc) {
  return(pnorm(a, ccc[[1]], ccc[[2]]))
}

fit_BP<-function(a){
  ccc <- as.vector(estimateBP(as.numeric(a), para.num = 4)$par)
  return(ccc)
}


### please be aware, I did some filtering on a when calling cdf ###
cdf_BP <- function(a, ccc, ngj=100) {
  a <- a[which(a > 0)]
  a <- sort(a)
  a <- a / ccc[4]
  a <- sort(x) / ccc[4]
  tt <- gauss.quad(ngj, kind = "jacobi", alpha = ccc[1] - 1, beta = ccc[2] - 1)
  pp <- rep(0, length(a))
  for (i in 1:length(a)) {
    pp[i] <- 1  /beta(ccc[1], ccc[2]) * 1 / 2 ^ (ccc[1] + ccc[2] - 1) * sum(ppois(a[i], (tt[[1]] + 1) * ccc[3] / 2 * (1 + tt[[1]]) ^ (ccc[1] - 1) * (1 - tt[[1]]) ^ (ccc[2] - 1)) * tt[[2]])
  }
  return(pp)
}

fit_P <- function(a) {
  ccc <- mean(a)
  return(ccc)
}

cdf_P <- function(a, ccc) {
  return(ppois(a, ccc))
}

#' @importFrom pscl zeroinfl
fit_ZIP <- function(x) {
  m1 <- pscl::zeroinfl(x ~ 1, dist = "poisson", EM = TRUE)
  pi_ML = predict(m1, type = "zero")[1]
  mean_ML = predict(m1, type = "count")[1]
  ccc <- c(pi_ML, mean_ML)
  return(ccc)
}

cdf_ZIP <- function(a, ccc) {
  pp <- ccc[1] + ppois(a, ccc[2]) * (1 - ccc[1])
  return(pp)
}

#' @importFrom MASS glm.nb
fit_NB <- function(x) {
  m1 <- MASS::glm.nb(x ~ 1)
  theta_ML = m1$theta
  mu_ML = exp(m1$coefficients[1])
  var_ML = mu_ML + (mu_ML ^ 2 / theta_ML)
  ccc <- c(theta_ML, mu_ML)
  return(ccc)
}

cdf_NB <- function(a, ccc) {
  pp <- pnbinom(a, mu = ccc[2], size = ccc[1])
  return(pp)
}

#' @importFrom pscl zeroinfl
fit_ZINB <- function(x) {
  m1 <- pscl::zeroinfl(x ~ 1, dist = "negbin", EM = TRUE)
  pi_ML = predict(m1, type = "zero")[1]
  theta_ML = m1$theta
  mean_ML = predict(m1, type = "count")[1]
  var_ML = mean_ML + (mean_ML ^ 2 / theta_ML)
  ccc <- c(pi_ML, theta_ML, mean_ML)
  return(ccc)
}

cdf_ZINB <- function(a, ccc) {
  pp <- ccc[1] + pnbinom(a, mu = ccc[3], size = ccc[2]) * (1 - ccc[1])
  return(pp)
}

#' @importFrom mclust Mclust
#' @importFrom mclust mclustBIC
fit_MG <- function(x) {
  aa <- mclust::Mclust(x, verbose = FALSE)
  p <- aa$parameters$pro
  m <- aa$parameters$mean
  s <- sqrt(aa$parameters$variance$sigmasq)
  if (length(s) != length(p)) {
    s <- rep(s, length(p))
  }
  ccc <- rbind(p, m, s)
  return(ccc)
}

cdf_MG <- function(a, ccc) {
  pp <- rep(0, length(a))
  for (i in 1:ncol(ccc)) {
    pp <- pp + ccc[1, i] * pnorm(a, ccc[2, i], ccc[3, i])
  }
  return(pp)
}

fit_ZIlogG <- function(x) {
  xx <- x
  x <- log(xx[which(xx != 0)])
  p0 <- mean(xx != 0)
  p <- 1 - p0
  m <- mean(x)
  s <- sd(x)
  ccc <- rbind(p, m, s)
  ccc <- cbind(c(p0, 0, 0), ccc)
  return(ccc)
}

cdf_ZIlogG <- function(a, ccc) {
  pp <- ccc[1, 1] + pnorm(log(a), ccc[2, 2], ccc[3, 2]) * ccc[1, 2]
  return(pp)
}

fit_ZIG <- function(x) {
  xx <- x
  x <- xx[which(xx != 0)]
  p0 <- mean(xx != 0)
  p <- 1 - p0
  m <- mean(x)
  s <- sd(x)
  ccc <- rbind(p, m, s)
  ccc <- cbind(c(p0, 0, 0), ccc)
  return(ccc)
}

cdf_ZIG <- function(a, ccc) {
  pp <- ccc[1, 1] + pnorm(a, ccc[2, 2], ccc[3, 2]) * ccc[1, 2]
  return(pp)
}

#' @importFrom mclust Mclust
#' @importFrom mclust mclustBIC
fit_ZIlogMG <- function(x) {
  xx <- x
  x <- log(xx[which(xx != 0)])
  p0 <- mean(xx != 0)
  aa <- mclust::Mclust(x, verbose = FALSE)
  p <- aa$parameters$pro
  p <- p * (1 - p0)
  m <- aa$parameters$mean
  s <- sqrt(aa$parameters$variance$sigmasq)
  if (length(s) != length(p)) {
    s <- rep(s, length(p))
  }
  ccc <- rbind(p, m, s)
  ccc <- cbind(c(p0, 0, 0), ccc)
  return(ccc)
}

cdf_ZIlogMG <- function(a, ccc) {
  pp <- rep(0, length(a))
  pp <- pp + ccc[1, 1]
  for (i in 2:ncol(ccc)) {
    pp <- pp + ccc[1, i] * pnorm(log(a), ccc[2, i], ccc[3, i])
  }
  return(pp)
}

#' @importFrom mclust Mclust
#' @importFrom mclust mclustBIC
fit_ZIMG <- function(x) {
  xx <- x
  x <- xx[which(xx != 0)]
  p0 <- mean(xx != 0)
  aa <- mclust::Mclust(x, verbose = FALSE)
  p <- aa$parameters$pro
  p <- p * (1 - p0)
  m <- aa$parameters$mean
  s <- sqrt(aa$parameters$variance$sigmasq)
  if (length(s) != length(p)) {
    s <- rep(s, length(p))
  }
  ccc <- rbind(p, m, s)
  ccc <- cbind(c(p0, 0, 0), ccc)
  return(list(result = ccc, k = aa$d))
}

cdf_ZIMG <- function(a, ccc) {
  pp <- rep(0, length(a))
  pp <- pp + ccc[1, 1]
  for (i in 2:ncol(ccc)) {
    pp <- pp + ccc[1, i] * pnorm(a, ccc[2, i], ccc[3, i])
  }
  return(pp)
}
