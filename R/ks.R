showmodel <- function(model, verbose = interactive()) {
  if (verbose) {
    cat(model)
    cat("\n")
  }
}

ks <- function(a, model, alternative = c("two.sided", "less", "greater"), verbose = FALSE) {
  showmodel(model, verbose)
  if (is.matrix(a)) {
    r <- ks.matrix(a, model, alternative)
  } else {
    r <- ks.vector(a, model, alternative)
  }
  return(r)
}

ks.matrix <- function(a, model, alternative = c("two.sided", "less", "greater"), verbose = FALSE) {
  apply(a, 1, function (x) ks.vector(x, model, alternative))
}

#' @importFrom kolmim pkolm
ks.vector <- function(a, model, alternative = c("two.sided", "less", "greater"), verbose = FALSE) {
  if (any(a == Inf)) return(list(ks = c("KS.stat" = Inf, "KS.pval" = 0), "fit_c" = NA, peak.count = -Inf))
  cdf_c = eval(parse(text=paste("cdf", model, sep = "_")))
  fit.result <- fit.distribution(a, model)
  if(is.null(fit.result) || is.null(fit.result$result)) return (list(ks = c("KS.stat" = Inf, "KS.pval" = 0), "fit_c" = NA, peak.count = -Inf))
  fit_c <- fit.result$result
  if (model %in% c("P", "NB", "G", "MG")) {
    r <- ks_statistics(a, fit_c, cdf_c)
  } else if (model %in% c("ZIP", "ZINB", "ZIG", "ZIlogG", "ZIMG", "ZIlogMG")) {
    r <- ks_statistics_ZI(a, fit_c, cdf_c)
  } else if (model %in% c("LTMG", "LTG")) {
    if (isTRUE(sum(abs(fit_c) == 0))) {
      return(list(ks = c("KS.stat" = 1, "KS.pval" = 0), "fit_c" = NA, peak.count = -Inf))
    }
    r <- ks_statistics_Zcut(a, fit_c, cdf_c)
  } else if (model %in% c("BP")) {
    warning("Something wrong with BP ")
    return(list(ks = c("KS.stat" = 1, "KS.pval" = 0), "fit_c" = NA, peak.count = -Inf))
  } else {
    browser()
  }
  n <- r$n
  x <- r$y - (0:(n - 1)) / n
  alternative <- match.arg(alternative)
  STATISTIC <- switch(alternative, two.sided = max(c(x, 1 / n - x)), greater = max(1 / n - x), less = max(x))
  ks_cut <- ks_determination(n)
  if (isTRUE(STATISTIC <  ks_cut)) {
    PVAL <- 1 - kolmim::pkolm(STATISTIC, n)
  } else {
    PVAL <- 0
  }
  return(list(ks = c("KS.stat" = STATISTIC, "KS.pval" = PVAL), "fit_c" = fit_c, peak.count = fit.result$peak.count))
}

ks_statistics <- function(a, fit_c, cdf_c = cdf_norm) {
  n <- length(a)
  x <- cdf_c(sort(a), fit_c)
  return (list("n" = n, "y" = x))
}

ks_statistics_ZI <- function(a, fit_c, cdf_c = cdf_ZIG) {
  aa1 <- a[which(a > 0)]
  p <- mean(a == 0)
  n <- length(aa1)
  x <- (cdf_c(sort(aa1), fit_c) - p) / (1 - p)
  return (list("n" = n, "y" = x))
}

ks_statistics_Zcut <- function(a, fit_c, cdf_c = cdf_LTMG1) {
  Zcut <- min(a[which(a != min(a))])
  aa <- a
  aa[which(aa < Zcut)] <- Zcut - 1
  aa1 <- aa[which(aa >= Zcut)]
  p <- mean(aa < Zcut)
  n <- length(aa1)
  x <- (cdf_c(sort(aa1), fit_c) - p) / (1 - p)
  return (list("n" = n, "y" = x))
}

#' @importFrom kolmim pkolm
ks_determination <- function(n) {
  pp <- min(1 / n, 0.1)
  pp_unit <- pp
  N <- 0
  ss <- 0
  while(ss < 0.5) {
    N <- N + 1
    STAT <- pp_unit * N
    ss <- kolmim::pkolm(STAT, n)
  }
  return(STAT)
}

#' M3S Fit
#'
#' @param x
#'
#' @param normalization
#' @param distribution models
#'
#' @export
M3Sfit <- function(x, normalization = c("none", "logplus", "log", "cpm.none", "cpm.logplus", "cpm.log"), distribution = c("G", "MG", "LTMG", "LTG")) {
  return (fit.distribution(normalize(x, normalization), distribution))
}

fit.distribution <- function(x, distribution = c("G", "MG", "LTMG", "LTG")) {
  if (any(x == Inf)) return()
  fitf_c = eval(parse(text=paste("fit", distribution, sep = "_")))
  if (distribution %in% c("P", "NB", "G", "MG", "ZIP", "ZINB", "ZIG", "ZIlogG", "ZIMG", "ZIlogMG")) {
    result <- fitf_c(x)
  } else if (distribution %in% c("LTMG", "LTG")) {
    Zcut <- min(x[which(x != min(x))])
    aa <- x
    aa[which(aa < Zcut)] <- Zcut - 1
    result <- fitf_c(aa, Zcut)
  } else {
    return()
  }
  if(!is.list(result)) {
    result <- list(result = result, peak.count = 1)
  }
  return (result)
}
