library(Rcpp)
library(RcppArmadillo)
library(inline)

#' Compute critical values for the (augmented) Dickey-Fuller test based on MacKinnon (2010).
#' @param n The sample size (numeric).
#' @param Type The deterministic part in the model ("none", "const", "trend").
#' @return A numeric vector with the 1%, 5%, 10% critical values.
#' @references MacKinnon (2010) Critical Values for Cointegration Tests, QED Working Paper 1227.
dfcv <- function(n, type = c("none", "const", "trend")) {
  if (!is.numeric(n) || n < 1) stop("n must be a positive number")
  type <- match.arg(type[1], c("none", "const", "trend"))
  switch(type,
         none =  c("1%" = -2.5658 - 1.960/n - 10.04/n^2,
                   "5%" = -1.9393 - 0.398/n,
                   "10%"= -1.6156 - 0.181/n),
         const = c("1%" = -3.4336 - 5.999/n - 29.25/n^2,
                   "5%" = -2.8621 - 2.738/n -  8.36/n^2,
                   "10%"= -2.5671 - 1.438/n -  4.48/n^2),
         trend = c("1%" = -3.9638 - 8.353/n - 47.44/n^2,
                   "5%" = -3.4126 - 4.039/n - 17.83/n^2,
                   "10%"= -3.1279 - 2.418/n -  7.58/n^2)
  )
}

#' Augmented Dickey-Fuller test
#' @param y The vector with the time series (numeric)
#' @param Type The deterministic part in the model ("none", "const", "trend")
#' @param lags The number of lags of the differenced time series
#' @return A list with three slots:
#' \itemize{
#' \item{$tau} the value of the statistic (tau)
#' \item{$criticals} the 1%, 5%, 10% critical values
#' \item{$regression} a matrix with three columns: coefficients, stderr, tstat
#' \item{$lm.fit.out} the output of lm.fit() used to run the regression
#' }
adf <- function(y, type = c("none", "const", "trend"), lags = 0) {
  if (!is.vector(y)) y <- as.numeric(y)
  type <- match.arg(type[1], c("none", "const", "trend"))
  n <- length(y)
  dy   <- y[-1] - y[-n]
  ylag <- y[-n]
  if (lags == 0) {
    switch(type,
           none  = {X <- matrix(ylag, n - 1, 1, dimnames = list(NULL, "y_1"))},
           const = {X <- cbind(y_1 = ylag, const = 1)},
           trend = {X <- cbind(y_1 = ylag, const = 1, trend = 1:(n-1))})
    ndx   <- complete.cases(cbind(dy, X))
    reg   <- lm.fit(X[ndx, , drop = FALSE], dy[ndx])
    rss   <- sum(reg$residuals^2)
    sig2  <- rss / reg$df.residual
    XXinv <- chol2inv(reg$qr$qr[drop = FALSE])
    se    <- sqrt(diag(XXinv) * sig2)
    structure(list(adf = c(tau = unname(reg$coefficients[1] / se[1])),
                   # criticals  = dfcv(n - 1 - lags, type),
                   lags = lags,
                   type = type,
                   regression = cbind(coeff = reg$coefficients, stderr = se, tstat = reg$coefficients/se),
                   lm.fit.out = c(reg, RSS = rss, sigma2 = sig2, list(vcov = XXinv * sig2))),
              class = "adf"
    )
  } else if (lags > 0) {
    switch(type,
           none  = {X <- cbind(ylag[-(1:lags)], embed(dy[-(n-1)], lags))
                    colnames(X) <- c("y_1", paste0("dy_", 1:lags))},
           const = {X <- cbind(ylag[-(1:lags)], embed(dy[-(n-1)], lags), const = 1)
                    colnames(X) <- c("y_1", paste0("dy_", 1:lags), "const")},
           trend = {X <- cbind(ylag[-(1:lags)], embed(dy[-(n-1)], lags), const = 1, trend = 1:(n - 1 - lags))
                    colnames(X) <- c("y_1", paste0("dy_", 1:lags), "const", "trend")}
           )
    dyshort <- dy[-(1:lags)]
    ndx   <- complete.cases(cbind(dyshort, X))
    reg   <- lm.fit(X[ndx, , drop = FALSE], dyshort[ndx])
    rss   <- sum(reg$residuals^2)
    sig2  <- rss / reg$df.residual
    XXinv <- chol2inv(reg$qr$qr[drop = FALSE])
    se    <- sqrt(diag(XXinv) * sig2)
    structure(list(adf = c(tau = unname(reg$coefficients[1] / se[1])),
                   criticals  = dfcv(n - 1 - lags, type),
                   lags = lags,
                   type = type,
                   regression = cbind(coeff = reg$coefficients, stderr = se, tstat = reg$coefficients/se),
                   lm.fit.out = c(reg, RSS = rss, sigma2 = sig2, list(vcov = XXinv * sig2))),
              class = "adf"
    )
  }
}

# logLik.adf <- function (object, REML = FALSE, ...) 
# {
#   res <- object$lm.fit.out$residuals
#   p   <- object$lm.fit.out$rank
#   N   <- length(res)
#   N0  <- N
#   if (REML) 
#     N <- N - p
#   val <- 0.5 * (- N * (log(2 * pi) + 1 - log(N) + log(sum(res^2))))
#   if (REML) 
#     val <- val - sum(log(abs(diag(object$lm.fit.out$qr$qr)[1L:p])))
#   attr(val, "nall") <- N0
#   attr(val, "nobs") <- N
#   attr(val, "df") <- p + 1
#   class(val) <- "logLik"
#   val
# }

#' Log-likelihood of the regression of the ADF test
#' @param object An object of class adf
#' @return An object of class \code{logLik} with \code{nobs} (number of observations),
#'         \code{df} (degrees of freedom) as attrubutes.
logLik.adf <- function(object) {
  n    <- length(object$lm.fit.out$residuals)
  RSS  <- object$lm.fit.out$RSS
  val <- -n/2 * (1 + log(2 * pi) - log(n) + log(RSS))
  attr(val, "nobs") <- n
  attr(val, "df") <- object$lm.fit.out$rank
  class(val) <- "logLik"
  val
}

#' Generic function that computes the correct Akaike Information Criterion.
#' The default function works if there exists a \code{logLik()} function that
#' returns a \code{logLik} object with attributes \code{df} and \code{nobs}.
#' @param object An object that can be passed to the \code{logLik()} function.
AICC <- function(object) UseMethod("AICC")

#' @describeIn AICC Default function for computing AICC.
AICC.default <- function(object) {
  ll <- logLik(object)
  p <- attr(ll, "df")
  n <- attr(ll, "nobs")
  2*(-as.numeric(ll) + p*n/(n - p - 1))
}

#' Generic function that computes the correct Hannan-Quinn Information Criterion.
#' The default function works if there exists a \code{logLik()} function that
#' returns a \code{logLik} object with attributes \code{df} and \code{nobs}.
#' @param object An object that can be passed to the \code{logLik()} function.
HQC <- function(object) UseMethod("HQC")

#' @describeIn HQC Default function for computing Hannan-Quinn Criterion.
HQC.default <- function(object) {
  ll <- logLik(object)
  p <- attr(ll, "df")
  n <- attr(ll, "nobs")
  2*(-as.numeric(ll) + p * log(log(n)))
}

#' Augmented Dickey-Fuller test with automatic lag selection.
#' @param y The vector with the time series (numeric).
#' @param Type The deterministic part in the model ("none", "const", "trend").
#' @param lags The maximum number of lags of the differenced time series.
#' @param selectlags The criterion for selecting the lags in the ADF regression.
#'        The choices are: "Fixed" (no lag selection), "AIC", "BIC", "AICC", "HQC".
#' @return An object of class \code{adftest}, which extends the class \code{adf},
#' with the slots:
#' \itemize{
#' \item{$tau} the value of the statistic (tau);
#' \item{$criticals} the 1%, 5%, 10% critical values;
#' \item{$regression} a matrix with three columns: coefficients, stderr, tstat;
#' \item{$lm.fit.out} the output of lm.fit() used to run the regression;
#' \item{$IC} a matrix with the information criterion values.
#' }
adf_test <- function(y, type = c("none", "const", "trend"), lags = 0,
                     selectlags = c("Fixed", "AIC", "BIC", "AICC", "HQC")) {
  if (!is.vector(y)) y <- as.numeric(y)
  type <- match.arg(type[1], c("none", "const", "trend"))
  selectlags <- match.arg(selectlags[1], c("Fixed", "AIC", "BIC", "AICC", "HQC"))
  
  if (selectlags == "Fixed" || lags == 0) return(adf(y, type, lags))

  outs      <- vector("list", lags + 1)
  aics      <- numeric(lags + 1)
  for (i in 0:lags) {
    if (i == lags) aics[i + 1] <- eval(parse(text = selectlags))(adf(y, type, i))
    else aics[i + 1] <- eval(parse(text = selectlags))(adf(y[-seq.int(lags-i)], type, i))
  }
  bestndx <- which.min(aics)
  rnames <- seq.int(0, lags)
  rnames[bestndx] <- paste0(rnames[bestndx], "*")
  adf <- adf(y, type, bestndx - 1)
  adf$IC <- matrix(aics, ncol = 1, dimnames = list(rnames, selectlags))
  class(adf) <- c("adftest", "adf")
  adf
}

#' Print an \code{adf} object.
#' @param x an object of \code{adf} class.
#' @return NULL
print.adf <- function(x) {
  n <- NROW(x$lm.fit.out$residuals)
  cat("Augmented Dickey-Fuller test\n")
  cat("Lag order (of AR representation) =", x$lags, "\n")
  cat("Length of time series =", n, "\n")
  cat(switch(x$type,
             none  = "No intercept or trend",
             const = "Intercept",
             trend = "Intercept and linear"
             ),
      "\n")
  cat("\nTau statistic =", x$adf, "\n")
  cat("Critical values:\n")
  print(dfcv(n - x$lm.fit.out$rank, x$type))
  cat("\nAugmente Dickey-Fuller regression:\n")
  print(cbind(x$regression, round(c(NA, 2*pnorm(-abs(x$regression[-1, 3]))), 4)))
}

#' Compute ADF on each column of y
#' @param y The vector with the time series (numeric).
#' @param Type The deterministic part in the model ("none", "const", "trend").
#' @param lags The maximum number of lags of the differenced time series.
#' @param selectlags The criterion for selecting the lags in the ADF regression.
#'        The choices are: "Fixed" (no lag selection), "AIC", "BIC", "AICC", "HQC".
adfs <- function(y, type = c("none", "const", "trend"), lags = 10, selectlags = "AIC") {
  # structure(sapply(as.data.frame(y), function(y) urca::ur.df(y, type, lags, selectlags)@teststat[1]),
  #           criticals = urca::ur.df(as.data.frame(x)[, 1], type, lags = 0)@cval[1, ])
  structure(sapply(as.data.frame(y), function(x) adf_test(x, type, lags, selectlags)$adf),
            criticals = dfcv(NROW(y), type))
}

#' Density of multivariate normal distribution.
#' @param x numeric vector at which the density is computed.
#' @param mean numeric vector of means.
#' @param cov numeric variance-covariance matrix.
#' @param log logical scalar, if \code{TRUE} log of density is computed, it defaults to \code{FALSE}.
#' @return A numeric scalar with the density or log-density value.
dmvn <- function(x, mean, cov, log = FALSE) {
  k <- length(mean)
  if (is.vector(x)) {
    if ((length(x) != k) || (ncol(cov) != k) || (nrow(cov) != k)) stop("check arguments' dimensions")
    x0 <- x - mean
  }
  else if (is.matrix(x)) {
    if ((ncol(x) != k) || (ncol(cov) != k) || (nrow(cov) != k)) stop("check arguments' dimensions")
    x0 <- t(x) - mean
  }
  else stop("x must be either a vector or a matrix")
  
  if (log) {
    dt <- determinant(cov, logarithm = TRUE)$modulus
    -0.5 * (k * log(2 * pi) + dt +  colSums(x0 * (chol2inv(chol(cov)) %*% x0)))
  } else {
    dt <- det(cov)
    (2 * pi)^(-0.5 * k) * dt^(-0.5) * exp(-0.5 * colSums(x0 * (chol2inv(chol(cov)) %*% x0)))
  }
}


#' Fit a vector autoregression.
#' @param y an object that can be converted into a matrix with the variables in the columns
#'        and observations in the rows.
#' @param type string that determines the deterministic part in the model: the choices are "none" (nothing),
#'        "const" (constant/intercept), "trend" (constant/intercept and linear trend).
#' @param lags the number of lags of the endogenous variables (VAR order).
#' @param x an object that can be converted into a matrix with exogenous variables
#'        in the columns.
#' @return An object of class \code{fitted_var}, which extends the class \code{var} with
#'         slots:
#' \itemize{
#' \item{coefficients} matrix of coefficients.
#' \item{residuals}	matrix of residuals.
#' \item{fitted.values}	matrix of fitted values.
#' \item{effects} matrix of orthogonal single-df effects.
#'                The first rank of them correspond to non-aliased coefficients,
#'                and are named accordingly.
#' \item{rank} integer, rank of design matrix.	
#' \item{df.residual} degrees of freedom of residuals.
#' \item{qr} the QR decomposition, see qr.
#' 
#' \item{cov.residuals} residual covariance matrix.
#' \item{lags} integer number of lags in the model (VAR order).
#' \item{nseries} integer number of time series.
#' \item{loglik} scalar value of log-likelihood.
#' \item{x} matrix of exogenous regressors used in the regressions.
#' \item{y} matrix of endogenous regressors used in the regressions.
#' }
var_fit <- function(y, type = c("none", "const", "trend"), lags = 1, x = NULL) {
  type = match.arg(type[1], c("none", "const", "trend"))
  y <- as.matrix(y)
  k <- ncol(y)
  n <- nrow(y)
  onek <- seq_len(k)
  onep <- seq_len(lags)
  if (!is.null(x)) {
    x <- as.matrix(x)
    if (n != nrow(x)) stop("y and x must have the same number of rows")
    X <- x[-onep, , drop = FALSE]
  } else X <- NULL
  if (is.null(k) || k < 2) stop("y must me a matrix-like object with more than one column")
  if (lags < 1) stop("lags must be a positive integer")
  switch(type,
         const = {X <- cbind(X, const = rep(1, n - lags))},
         trend = {X <- cbind(X, const = rep(1, n - lags), trend = seq.int(n - lags))}
  )
  Y <- embed(y, lags + 1)
  nms <- if (is.null(colnames(y))) paste0("y", onek) else colnames(y)
  colnames(Y) <- c(nms, paste(nms, rep(onep, each = k), sep="_"))
  keep <- complete.cases(Y, X)
  regx <- cbind(Y[keep, -onek, drop = F], X[keep, , drop = FALSE])
  regy <- Y[keep, onek]
  reg <- lm.fit(regx, regy)
  reg$cov.residuals <- crossprod(reg$residuals) / nrow(reg$residuals)
  reg$lags <- lags
  reg$nseries <- k
  reg$loglik <- sum(dmvn(reg$residuals, rep(0, k), reg$cov.residuals, TRUE))
  reg$x <- regx
  reg$y <- regy
  structure(reg, class = c("fitted_var", "var"))
}

#' Select the number of lags in a vector autoregression based on variuos criteria.
#' @param y an object that can be converted into a matrix with the variables in the columns
#'        and observations in the rows.
#' @param type string that determines the deterministic part in the model: the choices are "none" (nothing),
#'        "const" (constant/intercept), "trend" (constant/intercept and linear trend).
#' @param lags the number of lags of the endogenous variables (VAR order).
#' @param x an object that can be converted into a matrix with exogenous variables
#'        in the columns.
#' @return A list with two slots:
#' \itemize{
#' \item{IC} matrix with lags in the rows and criteria in the columns.
#' \item{best} named integer vector with lags selected by the various criteria.
#' }
var_lag_select <- function(y, type = c("none", "const", "trend"), lags = 4, x = NULL) {
  if (lags < 1) stop("pmax must be an integer greater than zero")
  y <- as.matrix(y)
  k <- ncol(y)
  n <- nrow(y) - lags
  x <- if (is.null(x)) NULL else as.matrix(x)
  ics <- matrix(NA_real_, lags, 6, dimnames = list(1:lags, c("loglik", "p(LR)", "AIC", "BIC", "HQC", "AICC")))
  for (i in seq.int(lags)) {
    if (i == lags) vfit <- var_fit(y, type, i, x)
    else {
      ndx <- seq.int(lags-i)
      vfit <- var_fit(y[-ndx, ], type, i, if (!is.null(x)) x[-ndx, ])
    }
    ll <- vfit$loglik
    ncoeff <- length(vfit$coefficients)
    ics[i, 1] <- ll
    if (i > 1) ics[i, 2] <- pchisq(2 * (ics[i, 1] - ics[i - 1, 1]), k^2, lower.tail = FALSE)
    ics[i, 3] <- 2*(-ll + ncoeff)/n
    ics[i, 4] <- (-2 * ll + log(n) * ncoeff)/n
    ics[i, 5] <- 2*(-ll + ncoeff * log(log(n)))/n
    ics[i, 6] <- 2*(-ll + ncoeff * n / (n - ncoeff - 1))/n
  }
  list(IC = ics, best = apply(ics[, -(1:2)], 2, which.min))
}

#' Compute Johansen's test for the selection of the cointegration rank.
#' @param y an object that can be converted into a matrix with the variables in the columns
#'        and observations in the rows.
#' @param type string that determines the deterministic part in the model: the choices are "none" (nothing),
#'        "rconst" (restricted constant), "uconst" (unrestricted constant),
#'        "rtrend" (restricted trend with unrestricted constant),
#'        "utrend" (unrestricted trend with unrestricted constant). Notice that *restricted* means that
#'        it enters only the cointegration relation(s).
#' @param lags the number of lags of the endogenous variables (VAR order).
#' @param rx an object that can be converted into matrix with the restricted regressors in the columns.
#' @param ux an object that can be converted into matrix with the unrestricted regressors in the columns.
#' @return An object of class \code{johansen} with slots
#' \itemize{
#' \item{statistics} a numeric matrix with eigenvalues, trace test and lmax test in the columns.
#' \item{nseries} integer with the number of endogenous time series in the system.
#' \item{lags} VAR order.
#' \item{n} integer with the length of each time series.
#' \item{type} type string that defines the deterministic part in the model (see the parameter above).
#' \item{nrx} integer with the number of restricted regressors (NULL if none).
#' \item{nux} integer with the number of unrestricted regressors (NULL if none).
#' }
#' list(statistics = out[kndx, ], nseries = k, lags = lags, n = n, type = type, nrx = ncol(rx), nux = ncol(ux))
johansen_test <- function(y, type = c("none", "rconst", "uconst", "rtrend", "utrend"), lags = 1,
                          rx = NULL, ux = NULL) {
  if (lags < 1) stop("lags must be an integer greater than zero")
  type <- match.arg(type[1], c("none", "rconst", "uconst", "rtrend", "utrend"))
  if (!is.matrix(y)) y <- as.matrix(y)
  if (!is.null(rx)) if (!is.matrix(rx)) rx <- as.matrix(rx)
  if (!is.null(ux)) if (!is.matrix(ux)) ux <- as.matrix(ux)
  n <- nrow(y)
  k <- ncol(y)
  dy <- y[-1, ] - y[-n, ]
  y_1 <- y[-n, ]
  if (type == "none" && lags == 1 && is.null(ux)) {
    keep <- complete.cases(dy, y_1)
    dyk <- dy[keep, , drop = FALSE]
    y_1k <- y_1[keep, , drop = FALSE]
    Sxx <- crossprod(y_1k, y_1k)
    Sxy <- crossprod(y_1k, dyk)
    Syy <- crossprod(dyk, dyk)
    lam <- geigen::geigen(Sxy %*% chol2inv(chol(Syy)) %*% t(Sxy), Sxx, symmetric = T, only.values = T)$values
    trace_stat <- -(n - 1) * cumsum(log(1 - lam))
    lmax_stat  <- -(n - 1) * log(1 - lam)
    out <- cbind(eigenvalue = lam, trace = trace_stat, lmax = lmax_stat)[length(lam):1, ]
    rownames(out) <- 0:(length(lam) - 1)
  } else {
    kndx <- seq.int(k)
    lndx <- seq.int(lags - 1)
    DY <- embed(dy, lags)
    RX <- cbind(y_1    = if (lags > 1) y_1[-lndx, ] else y_1,
                rconst = if (type == "rconst") rep(1, n - lags) else NULL,
                rtrend = if (type == "rtrend") seq.int(n - lags) else NULL,
                rx[-seq.int(lags), , drop = FALSE])
    UX <- cbind(dy_    = if (lags > 1) DY[, -kndx] else NULL,
                uconst = if (type == "uconst" || type == "rtrend" || type == "utrend") rep(1, n - lags) else NULL,
                utrend = if (type == "utrend") seq.int(n - lags) else NULL,
                ux[-seq.int(lags), , drop = FALSE])
    if (is.null(UX)) {
      keep <- complete.cases(DY, RX)
      DYk <- DY[keep, , drop = FALSE]
      RXk <- RX[keep, , drop = FALSE]
      Sxx <- crossprod(RXk, RXk)
      Sxy <- crossprod(RXk, DYk[, kndx])
      Syy <- crossprod(DYk[, kndx], DYk[, kndx])
      lam <- geigen::geigen(Sxy %*% chol2inv(chol(Syy)) %*% t(Sxy), Sxx,
                            symmetric = T, only.values = T)$values
      trace_stat <- -(n - lags) * cumsum(log(1 - lam))
      lmax_stat  <- -(n - lags) * log(1 - lam)
      out <- cbind(eigenvalue = lam, trace = trace_stat, lmax = lmax_stat)[length(lam):1, ]
      rownames(out) <- 0:(length(lam) - 1)
    } else {
      keep <- complete.cases(DY, RX, UX)
      DYres <- lm.fit(UX[keep, , drop = FALSE], DY[keep, kndx, drop = FALSE])$residuals
      RXres <- lm.fit(UX[keep, , drop = FALSE], RX[keep, , drop = FALSE])$residuals
      Sxx <- crossprod(RXres, RXres)
      Sxy <- crossprod(RXres, DYres)
      Syy <- crossprod(DYres, DYres)
      lam <- geigen::geigen(Sxy %*% chol2inv(chol(Syy)) %*% t(Sxy), Sxx,
                            symmetric = T, only.values = T)$values
      trace_stat <- -(n - lags) * cumsum(log(1 - lam))
      lmax_stat  <- -(n - lags) * log(1 - lam)
      out <- cbind(eigenvalue = lam, trace = trace_stat, lmax = lmax_stat)[length(lam):1, ]
      rownames(out) <- 0:(length(lam) - 1)
    }
  }
  structure(list(statistics = out[kndx, ], nseries = k, lags = lags, n = n, type = type, nrx = ncol(rx), nux = ncol(ux)),
            class = "johansen")
}

#' \code{print} method for the \code{johansen} class.
#' @param x an object of \code{johansen} class.
#' @return NULL
print.johansen <- function(x) {
  cat("Johansen test\n")
  cat("Number of equations =", x$nseries, "\n")
  cat("Lag order (of VAR representation) =", x$lags, "\n")
  cat("Length of time series =", x$n, "\n")
  cat(switch(x$type,
             none = "Case 1: No constant or trend",
             rconst = "Case 2: Constant in cointegration vector",
             uconst = "Case 3: Unrestricted constant",
             rtrend = "Case 4: Unrestricted constant and trend in cointegration vector",
             utrend = "Case 5: Unrestricted constant and trend"),
      "\n")
  
  cv <- switch(x$type,
               none = NULL,
               rconst = list(lmax = matrix(c(7.52, 13.75, 19.77, 25.56, 31.66, 37.45, 43.25, 48.91, 54.35, 60.25,
                                             66.02, 9.24, 15.67, 22, 28.14, 34.4, 40.3, 46.45, 52, 57.42, 63.57,
                                             69.74, 12.97, 20.2, 26.81, 33.24, 39.79, 46.82, 51.91, 57.95, 63.71,
                                             69.94, 76.63),
                                           11, 3, dimnames = list(1:11, c("10%", "5%", "1%"))),
                             trace = matrix(c(7.52, 17.85, 32, 49.65, 71.86, 97.18, 125.58, 159.48, 196.37, 236.54,
                                              282.45, 9.24, 19.96, 34.91, 53.12, 76.07, 102.14, 131.7, 165.58,
                                              202.92, 244.15, 291.4, 12.97, 24.6, 41.07, 60.16, 84.45, 111.01,
                                              143.09, 177.2, 215.74, 257.68, 307.64),
                                            11, 3, dimnames = list(1:11, c("10%", "5%", "1%")))
               ),
               uconst = list(lmax = matrix(c(2.69, 12.07, 18.6, 24.73, 30.9, 36.76, 42.32, 48.33, 53.98, 59.62,
                                             65.38, 3.76, 14.07, 20.97, 27.07, 33.46, 39.37, 45.28, 51.42, 57.12,
                                             62.81, 68.83, 6.65, 18.63, 25.52, 32.24, 38.77, 45.1, 51.57, 57.69,
                                             62.8, 69.09, 75.95),
                                           11, 3, dimnames = list(1:11, c("10%", "5%", "1%"))),
                             trace = matrix(c(2.69, 13.33, 26.79, 43.95, 64.84, 89.48, 118.5, 150.53, 186.39,
                                              225.85, 269.96, 3.76, 15.41, 29.68, 47.21, 68.52, 94.15, 124.24, 156,
                                              192.89, 233.13, 277.71, 6.65, 20.04, 35.65, 54.46, 76.07, 103.18, 133.57,
                                              168.36, 204.95, 247.18, 293.44),
                                            11, 3, dimnames = list(1:11, c("10%", "5%", "1%")))
               ),
               rtrend = list(lmax = matrix(c(10.49, 16.85, 23.11, 29.12, 34.75, 40.91, 46.32, 52.16, 57.87, 63.18,
                                             69.26, 12.25, 18.96, 25.54, 31.46, 37.52, 43.97, 49.42, 55.5, 61.29,
                                             66.23, 72.72, 16.26, 23.65, 30.34, 36.65, 42.36, 49.51, 54.71, 62.46,
                                             67.88, 73.73, 79.23),
                                           11, 3, dimnames = list(1:11, c("10%", "5%", "1%"))),
                             trace = matrix(c(10.49, 22.76, 39.06, 59.14, 83.2, 110.42, 141.01, 176.67, 215.17,
                                              256.72, 303.13, 12.25, 25.32, 42.44, 62.99, 87.31, 114.9, 146.76,
                                              182.82, 222.21, 263.42, 310.81, 16.26, 30.45, 48.45, 70.05, 96.58,
                                              124.75, 158.49, 196.08, 234.41, 279.07, 327.45),
                                            11, 3, dimnames = list(1:11, c("10%", "5%", "1%")))
               ),
               utrend = list(lmax = matrix(c(2.57, 14.84, 21.53, 27.76, 33.74, 39.5, 45.49, 51.14, 57.01, 62.69,
                                             68.22, 3.74, 16.87, 23.78, 30.33, 36.41, 42.48, 48.45, 54.25, 60.29,
                                             66.1, 71.68, 6.4, 21.47, 28.83, 35.68, 41.58, 48.17, 54.48, 60.81,
                                             66.91, 72.96, 78.51),
                                           11, 3, dimnames = list(1:11, c("10%", "5%", "1%"))),
                             trace = matrix(c(2.57, 16.06, 31.42, 50.74, 73.4, 100.14, 130.84, 164.34, 201.95,
                                              244.12, 288.08, 3.74, 18.17, 34.55, 54.64, 77.74, 104.94, 136.61, 170.8,
                                              208.97, 250.84, 295.99, 6.4, 23.46, 40.49, 61.24, 85.78, 114.36, 146.99,
                                              182.51, 222.46, 263.94, 312.58),
                                            11, 3, dimnames = list(1:11, c("10%", "5%", "1%")))
               )
  )
  if (x$type ==  "none" || x$nseries > 11) {
    print(x$test)
  } else {
    cat("\nTrace test:\n")
    print(cbind(round(x$statistics[, 2, drop = FALSE], 2), cv$trace[x$nseries:1, ]))
    cat("\nLamx test:\n")
    print(cbind(round(x$statistics[, 3, drop = FALSE], 2), cv$lmax[x$nseries:1, ]))
  }
}

vecm_fit <- function(y, type = c("none", "rconst", "uconst", "rtrend", "utrend"), lags = 1, crank = 1,
                     rx = NULL, ux = NULL) {
  if (lags < 1) stop("lags must be an integer greater than zero")
  type <- match.arg(type[1], c("none", "rconst", "uconst", "rtrend", "utrend"))
  if (!is.matrix(y)) y <- as.matrix(y)
  if (!is.null(rx)) if (!is.matrix(rx)) rx <- as.matrix(rx)
  if (!is.null(ux)) if (!is.matrix(ux)) ux <- as.matrix(ux)
  n <- nrow(y)
  k <- ncol(y)
  if (crank < 1 || crank >= k) stop("cointegration rank must be a positive integer smaller than the number of endogenous time series")
  dy <- y[-1, ] - y[-n, ]
  colnames(dy) <- paste0("d_", colnames(y))
  y_1 <- y[-n, ]
  kndx <- seq.int(k)
  lndx <- seq.int(lags - 1)
  rndx <- seq.int(crank)
  if (type == "none" && lags == 1 && is.null(ux)) {
    Sxx <- crossprod(y_1, y_1)
    Sxy <- crossprod(y_1, dy)
    Syy <- crossprod(dy, dy)
    eig <- geigen::geigen(Sxy %*% chol2inv(chol(Syy)) %*% t(Sxy), Sxx, symmetric = T)
    B <- eig$vectors[, k:(k - crank + 1), drop = FALSE]
    B <- B %*% solve(B[rndx, rndx])
    colnames(B) <- paste0("EC", rndx)
    ec <- y_1 %*% B
    colnames(ec) <- paste0("EC", rndx)
    keep <- complete.cases(ec, dy)
    reg <- lm.fit(ec, dy[keep, ])
  } else {
    DY <- embed(dy, lags)
    colnames(DY) <- c(paste0("d_", colnames(y)),
                      if (lags > 1) paste0("d_", colnames(y), "_", rep(1:(lags - 1), each = k)))
    RX <- cbind(y_1    = if (lags > 1) y_1[-lndx, ] else y_1,
                rconst = if (type == "rconst") rep(1, n - lags) else NULL,
                rtrend = if (type == "rtrend") seq.int(n - lags) else NULL,
                rx[-seq.int(lags), , drop = FALSE])
    UX <- cbind(dy_    = if (lags > 1) DY[, -kndx] else NULL,
                uconst = if (type == "uconst" || type == "rtrend" || type == "utrend") rep(1, n - lags) else NULL,
                utrend = if (type == "utrend") seq.int(n - lags) else NULL,
                ux[-seq.int(lags), , drop = FALSE])
    keep <- complete.cases(DY, RX, UX)
    if (is.null(UX)) {
      Sxx <- crossprod(RX[keep, , drop = FALSE], RX[keep, , drop = FALSE])
      Sxy <- crossprod(RX[keep, , drop = FALSE], DY[, kndx, drop = FALSE])
      Syy <- crossprod(DY[keep, kndx, drop = FALSE], DY[keep, kndx, drop = FALSE])
      eig <- geigen::geigen(Sxy %*% chol2inv(chol(Syy)) %*% t(Sxy), Sxx, symmetric = T)
      B <- eig$vectors[, k:(k - crank + 1), drop = FALSE]
      B <- B %*% solve(B[rndx, rndx])
      colnames(B) <- paste0("EC", rndx)
      ec <- RX[keep, , drop = FALSE] %*% B
      colnames(ec) <- paste0("EC", rndx)
      reg <- lm.fit(ec, dy[keep, ])
    } else {
      DYres <- lm.fit(UX[keep, , drop = FALSE], DY[keep, kndx, drop = FALSE])$residuals
      RXres <- lm.fit(UX[keep, , drop = FALSE], RX[keep, , drop = FALSE])$residuals
      Sxx <- crossprod(RXres, RXres)
      Sxy <- crossprod(RXres, DYres)
      Syy <- crossprod(DYres, DYres)
      eig <- geigen::geigen(Sxy %*% chol2inv(chol(Syy)) %*% t(Sxy), Sxx, symmetric = T)
      B <- eig$vectors[, length(eig$values):(length(eig$values) - crank + 1), drop = FALSE]
      B <- B %*% solve(B[rndx, rndx])
      colnames(B) <- paste0("EC", rndx)
      ec <- RX[keep, , drop = FALSE] %*% B
      colnames(ec) <- paste0("EC", rndx)
      reg <- lm.fit(cbind(UX[keep, , drop = FALSE], ec), DY[keep, kndx, drop = FALSE])
    }
  }
  reg <- unclass(reg)
  reg$B <- B
  reg$EC <- ec
  reg$cov.residuals <- crossprod(reg$residuals) / nrow(reg$residuals)
  reg$nseries = k
  reg$lags = lags
  reg$crank = crank
  reg$type = type
  reg$nrx = ncol(rx)
  reg$nux = ncol(ux)
  reg$loglik <- sum(dmvn(reg$residuals, rep(0, k), reg$cov.residuals, TRUE))
  reg$y <- y
  structure(reg, class = c("fitted_vecm", "vecm"))
}

vecm2var <- function(vecm) {
  if (!("vecm" %in% class(vecm))) stop("argument must be of 'vecm' class")
  crank <- ncol(vecm$B)
  lags <- as.integer(vecm$lags)
  k    <- ncol(vecm$residuals)
  cf0  <- vecm$coefficients
  ncf0 <- nrow(cf0)
  type <- vecm$type
  B    <- vecm$B
  vnam <- rownames(B)[1:k]
  PIt  <- B %*% cf0[(ncf0 - crank + 1):ncf0, , drop = FALSE]
  nrx  <- nrow(B) - k               # number of restricted regressors
  nux  <- ncf0 - ncol(B) - k * (lags - 1) # number of unrestricted regressors
  
  if (lags == 1) {
    cf1 <- diag(k) + PIt[1:k, ]
    colnames(cf1) <- vnam
    rownames(cf1) <- paste0(vnam, "_", 1)
    if (nrx > 0) cf1 <- rbind(cf1, PIt[(k+1):nrow(PIt), , drop = FALSE]) # if there are restricted regressors
    if (nux > 0) cf1 <- rbind(cf1, cf0[1:nux, , drop = FALSE]) # if there are unrestricted regressors
  } else {
    cf1 <- matrix(0, k * lags, k, dimnames = list(paste0(vnam, "_", rep(1:lags, each = k)), vnam))
    cf1[1:k, ] <- diag(k) + PIt[1:k, ] + cf0[1:k, ]
    for (i in 2:lags) {
      if (i == lags) cf1[((i - 1) * k + 1):(i * k), ] <- -cf0[((i - 2) * k + 1):((i - 1) * k) , ]
      else cf1[((i - 1) * k + 1):(i * k), ] <- cf0[((i - 1) * k + 1):(i * k)] - cf0[((i - 2) * k + 1):((i - 1) * k) , ]
      if (nrx > 0) cf1 <- rbind(cf1, PIt[(k+1):nrow(PIt), , drop = FALSE]) # if there are restricted regressors
      if (nux > 0) cf1 <- rbind(cf1, cf0[((lags - 1)*k + 1):((lags - 1)*k + nux), , drop = FALSE]) # if there are unrestricted regressors
    }
  }
  colnames(vecm$cov.residuals) <- rownames(vecm$cov.residuals) <- vnam
  structure(list(coefficients = cf1, cov.residuals = vecm$cov.residuals,
                 y = vecm$y, rx = vecm$rx, ux = vecm$ux, lags = lags, nseries = k),
            class = "var")
}

var_sim <- function(var, err, init = 0) {
  if (!("var" %in% class(var))) stop("var must be of class 'var'")
  if (!is.matrix(err)) err <- as.matrix(err)
  if (ncol(err) != var$nseries) stop(paste("'err' must have", var$nseries, "columns"))
  n <- nrow(err)
  k <- var$nseries
  onek <- 1:k
  lags <- var$lags
  y <- matrix(0, n + lags, k, dimnames = list(-lags:(n - 1), colnames(var$coefficients)))
  y[1:lags, ] <- init
  phi <- var$coefficients[1:(k * lags), , drop = FALSE]
  for (t in seq.int(lags + 1, lags + n)) {
    y[t, ] <- embed(y[(t - lags):(t-1), , drop = FALSE], lags) %*% phi + err[t - lags, ]
  }
  y
}


# ---- C++ function to simulate pass a VAR(p) filter
# varFilter(mAt, mEt, mIt)
# INPUT:
# mAt : k*p x p matrix with the VAR(p) coefficients (traspoed)
# mEt : n x k matrix of time series to be filtered, observations in rows and variables in columns
# mIt : k*p x p matrix with initial values for the VAR(p) filter
# OUTPUT:
# It returns (p + n) x k matrix computed as follows:
# Yt[1:p, ] = mIt
# Yt[t, ] = Y[t-1, ] %*% mAt[1:k, ] + ... + Y[((p-1)*k + 1):p*k, ] %*% mAt[] + mEt[t, ]
sourceCpp(code='
          // [[Rcpp::depends(RcppArmadillo)]]
          
          #include <RcppArmadillo.h>
          
          // [[Rcpp::export]]
          arma::mat varFilter(arma::mat mAt, arma::mat mEt, arma::mat mIt) {
          int n = mEt.n_rows;
          int k = mEt.n_cols;
          int p = mAt.n_rows / k;
          arma::mat mYt;
          mYt.zeros(n + p, k);
          mYt.rows(0, p - 1) = mIt;
          for (int t = p; t < n + p; t++) {
            for (int i = 0; i < p; i++) {
              mYt.row(t) = mYt.row(t) + mYt.row(t - 1 - i) * mAt.rows(i*k, (i + 1)*k - 1);
            }
            mYt.row(t) = mYt.row(t) + mEt.row(t - p);
          }
          return mYt;
          }'
)


vecm_irf_boot <- function(fitted_vecm, horizon = 20, nboot = 100) {
  if (!("fitted_vecm" %in% class(fitted_vecm))) stop("first argument must be of class 'fitted_vecm'")
  var   <- vecm2var(fitted_vecm)
  init  <- fitted_vecm$y[complete.cases(fitted_vecm$y), ][1:fitted_vecm$lags, , drop = FALSE]
  res   <- fitted_vecm$residuals
  n     <- nrow(res)
  k     <- fitted_vecm$nseries
  lags  <- fitted_vecm$lags
  type  <- fitted_vecm$type
  crank <- fitted_vecm$crank
  rx    <- fitted_vecm$rx
  ux    <- fitted_vecm$ux
  imp   <- matrix(0, horizon, k)
  irf   <- vector("list", k)
  names(irf) <- colnames(fitted_vecm$y)
  Pt    <- chol(fitted_vecm$cov.residuals)
  zeros <- matrix(0, lags, k)
  cosel <- 1:(lags*k)
  for (i in 1:k) {
    imp[1, ] <- Pt[i, ]
#    irf[[i]] <-var_sim(vecm2var(fitted_vecm), imp)[-(1:lags), ]
    varred <- vecm2var(fitted_vecm)
    irf[[i]] <- varFilter(varred$coefficients[cosel, ], imp, zeros)[-(1:lags), ]
    colnames(irf[[i]]) <- names(irf)
    rownames(irf[[i]]) <- 0:(horizon - 1)
  }
  irfs  <- vector("list", k)
  names(irfs) <- colnames(fitted_vecm$y)
  for (i in 1:k) irfs[[i]] <- array(0, c(horizon, k, nboot),
                                    dimnames = list(0:(horizon - 1), names(irfs), NULL))
  for (i in 1:nboot) {
    eps <- res[sample.int(n, replace = TRUE), ]
#    y <- var_sim(var, eps, init)
    y <- varFilter(var$coefficients[cosel, ], eps, init)
    colnames(y) <- names(irfs)
    new_vecm <- vecm_fit(y, type, lags, crank, rx, ux)
    new_var <- vecm2var(new_vecm)
    Pt <- chol(new_vecm$cov.residuals)
    for (j in 1:k) {
      imp[1, ] <- Pt[j, ]
#      irfs[[j]][, , i] <- var_sim(var, imp)[-(1:lags), ]
      irfs[[j]][, , i] <- varFilter(new_var$coefficients, imp, zeros)[-(1:lags), ]
    }
  }
  structure(list(irf = irf, bootirf = irfs), class = "bootirf")
}

plot.bootirf <- function(x, probs = seq.int(5, 95, 10)/100, meanline = TRUE, ...) {
  probs <- as.numeric(probs)
  np <- length(probs)
  if (np/2 != np %/% 2) stop("second argument 'probs' must have an even number of elements")
  np2 <- np/2
  k <- length(x$irf)
  horizon <- nrow(x$irf[[1]])
  xgrid <- 0:(horizon - 1)
  nms <- colnames(x$irf[[1]])
  for (i in 1:k) {
    for (j in 1:k) {
      quans <- apply(x$bootirf[[i]][, j, ], 1, quantile, probs = probs)
      rng <- range(c(quans, 0))
      plot(x = xgrid, y = rep(0, horizon), xlab = "", ylab = "",
           ylim = rng, lty = 2, type = "l",
           main = paste(nms[i], "->", nms[j]), ...)
      for (v in 1:np2) {
        polygon(x = c(xgrid, rev(xgrid)), y = c(quans[v, ], rev(quans[np + 1 - v, ])),
                col = paste0("gray", 100 - v * 5), border = FALSE)
      }
      if (meanline) lines(x = xgrid, y = x$irf[[i]][, j], lwd = 2)
      readline(prompt="Press [enter] to continue")
    }
  }
}


# Function for extracting selection rates of Johansen's test
Johansen_select <- function(stat, crit) {
  k <- sum(grepl("k_star",colnames(stat)))
  n <- nrow(stat)
  less <- stat < rep(crit, each = n)
  val <- apply(cbind(less, TRUE), 1, function(x) {which(x)[1] - 1})
  df <- data.frame(table(factor(val, as.character(0:k))) / n * 100)
  colnames(df) <- c("Rels","Sel_rate")
  df <- df %>%
    mutate(Rels = paste0("r = ",0:k))
  df
}