is.discrete <- function(x) {
  return(!isFALSE(all(x %% 1 == 0)))
}

M3SObject <- setClass("M3SObject", slots = c(all.nonnegative = "logical", neg.inf = "logical", significant.zeros = "logical", is.discrete = "logical"))

M3SModelBase <- setClass("M3SModelBase", representation("VIRTUAL"))
setGeneric("fit", function(x) standardGeneric("fit"))
setGeneric("pdf", function(x) standardGeneric("pdf"))

M3SModelP <- setClass("M3SModelP", contains = "M3SModelBase")
setMethod("fit", "M3SModelBase", function(x) fit_P(x))
setMethod("pdf", "M3SModelBase", function(x) cdf_P(x))

CPM <- function(x) {
  if(is.matrix(x)) {
    r <- CPM.matrix(x)
  } else {
    r <- CPM.vector(x)
  }
  r[which(is.na(r))] <- 0
  return(r)
}

#' @importFrom plyr aaply
CPM.matrix <- function(x) {
  plyr::aaply(x, 1, CPM.vector)
}

CPM.vector <- function(x) {
  x * 10 ^ 6 / sum(x)
}

#' M3S
#'
#' @param x A matrix or an array stores a the data, if a matrix, each row is one feature and each column is one gene
#' @param sampling Unused for now
#' @param significant.cut Cut off for significant.zeros
#' @param verbose display promote or not
#'
#' @return statistical indicator for available modes
#' @export
M3S <- function(x, sampling = 100, significant.cut = 0.1, verbose = interactive()) {
  if(!is.matrix(x)) stop("Only matrix supported.")
  o <- M3SObject()
  # check if the data are all non-negative
  o@all.nonnegative <- all(x >= 0)
  # check if the data are with a significant number of 0s
  o@significant.zeros <- mean(x == 0) > significant.cut
  # if the data are discrete or continuous data
  o@is.discrete <- is.discrete(x)
  o@neg.inf <- any(x == -Inf)
  if (verbose) {
    print(o)
  }
  # If the number of samples is smaller than 30, we will no conduct the analysis and return sample size too small
  if (length(x) < 30) {
    warning("sample size too small")
    return(NA)
  }

  X_selected <- x
  M <- nrow(x)
  # If nrow of the data is larger than 10000, with out a specific setting by the user, we will do sampling of the data
  if (is.matrix(x) && M > 10000) {
    # Sampling method (1),
    cc <- apply(x, 1, mean)
    Selected_IDs <- c()
    for (K in 1:100) {
      Selected_IDs <- c(Selected_IDs, sample(floor(M / 100) * (K - 1) + 1:floor(M / 100) * K, 1))
    }
    X_selected <- x[order(cc)[Selected_IDs], ]
  }
  # Data fitting fo the X_selecte
  if (o@is.discrete && !o@all.nonnegative) {
    return(NA)
  } else if (o@is.discrete && o@significant.zeros) {
    models.none <- c("P", "BP", "NB", "ZIP", "ZINB")
    models.logplus <- c()
    models.log <- c()
    models.cpm.none <- c("G", "MG", "ZIG", "ZIlogG", "ZIlogMG", "ZIMG")
    models.cpm.logplus <- c("G", "MG", "ZIG", "ZIlogG", "ZIlogMG", "ZIMG")
    models.cpm.log <- c("LTMG", "LTG")
  } else if (o@is.discrete && !o@significant.zeros) {
    models.none <- c("P", "BP", "NB")
    models.logplus <- c()
    models.log <- c()
    models.cpm.none <- c("G", "MG")
    models.cpm.logplus <- c("G", "MG")
    models.cpm.log <- c("LTMG", "LTG")
  } else if (!o@is.discrete && !o@all.nonnegative && (!o@neg.inf)) {
    models.none <- c("G", "MG", "LTMG", "LTG")
    models.logplus <- c()
    models.log <- c()
    models.cpm.none <- c()
    models.cpm.logplus <- c()
    models.cpm.log <- c()
  } else if (!o@is.discrete && !o@all.nonnegative && (o@neg.inf)) {
    models.none <- c("LTMG", "LTG")
    models.logplus <- c()
    models.log <- c()
    models.cpm.none <- c()
    models.cpm.logplus <- c()
    models.cpm.log <- c()
  } else if (!o@is.discrete && o@all.nonnegative && o@significant.zeros) {
    models.none <- c("G", "MG", "LTMG", "LTG")
    models.logplus <- c("G", "MG", "ZIG", "ZIlogG", "ZIlogMG", "ZIMG")
    models.log <- c("LTMG", "LTG")
    models.cpm.none <- c("G", "MG", "LTMG", "LTG")
    models.cpm.logplus <- c("G", "MG", "ZIG", "ZIlogG", "ZIlogMG", "ZIMG")
    models.cpm.log <- c("LTMG", "LTG")
  } else if (!o@is.discrete && o@all.nonnegative && !o@significant.zeros) {
    models.none <- c("G", "MG")
    models.logplus <- c("G", "MG")
    models.log <- c()
    models.cpm.none <- c("G", "MG")
    models.cpm.logplus <- c("G", "MG")
    models.cpm.log <- c()
  } else {
    return(NA)
  }
  ret <- list()
  for (normalization in c("none", "logplus", "log", "cpm.none", "cpm.logplus", "cpm.log")) {
    print(paste("normalization", normalization, sep = ": "))
    for (model in eval(parse(text = paste0("models.", normalization)))) {
      ret[[paste(model, normalization, sep = ".")]] <- ks(normalize(X_selected, normalization), model, verbose = verbose)
    }
  }

  nret <- matrix(nrow = nrow(X_selected), ncol = length(ret), dimnames = list(rownames(X_selected), names(ret)))
  p <- matrix(nrow = nrow(X_selected), ncol = length(ret), dimnames = list(rownames(X_selected), names(ret)))

  peak.counts <- rep(0, length(ret))
  for (col in 1:length(ret)) {
    for (row in 1:nrow(X_selected)) {
      nret[row, col] <- ret[[col]][[row]]$ks[2]
      if(ret[[col]][[row]]$peak.count == -Inf){
        peak.counts[col] <- peak.counts[col] + 10
      } else {
        peak.counts[col] <- peak.counts[col] + ret[[col]][[row]]$peak.count
      }
    }
  }

  for (col in 1:length(ret)) {
    p[, col] <- p.adjust(nret[, col], method = "fdr")
  }

  candidate <- rep(NA, length(ret))

  length.less.than.0.1 <- rep(NA, length(ret))

  for (col in 1:length(ret)) {
    length.less.than.0.1[col] <- length(which(p[, col] < 0.1))
    candidate[col] <- length.less.than.0.1[col]< nrow(X_selected) / 10
  }

  if(all(!candidate)) {
    candidate[which.max(length.less.than.0.1)] <- TRUE
  }

  weights <- list(P = 1, G = 2, NB = 2.5, ZIP = 3, ZINB = 4.8, ZIG = 4.5, LTG = 4, BP = 5, MG = 6, ZIMG = 10, LTMG = 10)

  model.names <- unlist(lapply(strsplit(names(ret), "[.]"), `[[`, 1))
  selected.weight <- weights[model.names][candidate]
  minmodel <- min(unlist(selected.weight))
  minindex <- which.min(unlist(selected.weight))

  choosed.index <- 10

  best <- names(ret)[candidate][minindex]

  if(selected.weight[minindex] == choosed.index) {
    which.weight <- which(selected.weight == choosed.index)
    if(length(which.weight) > 1) {
      if (verbose) {
        cat("We must decide ZIMG or LTMG, :\n")
        print(peak.counts[candidate][which.weight])
        print(names(ret)[candidate][which.weight])
      }
      best <- names(ret)[candidate][which.weight][which.min(peak.counts[candidate][which.weight])]
    }
  }
  return(list(best.model = best, p = p, ret))
}

normalize <- function(x, normalization = c("none", "logplus", "log", "cpm.none", "cpm.logplus", "cpm.log")) {
  eval(parse(text = paste0("normalize.", normalization, "(x)")))
}

normalize.none <- function(x) {
  x
}

normalize.logplus <- function(x) {
  log(x + 1)
}

normalize.log <- function(x) {
  log(x) + 1
}

normalize.cpm.none <- function(x) {
  CPM(x)
}

normalize.cpm.logplus <- function(x) {
  log(CPM(x) + 1)
}

normalize.cpm.log <- function(x) {
  log(CPM(x)) + 1
}
