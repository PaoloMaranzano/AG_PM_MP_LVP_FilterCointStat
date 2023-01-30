# Extract the RW component using both the filter and the smoother from a RW + WN model.
# It returns a matrix with filtered (FLT) and smoothed (SMO) RW estimates.
trendextractor <- function(x) {
  n <- NROW(x)
  lvard <- log(var(x[-1] - x[-n])) # variance of differenced series
  mod <- SSModel(x ~ SSMtrend(1, NA), H = NA) # local level model
  obj <- function(pars) { # computes -meanLogLik at parameter values in pars
    mod$Q[1, 1, 1] <<- exp(pars[1])
    mod$H[1, 1, 1] <<- exp(pars[2])
    -logLik(mod) / n
  }
  opt <- optim(c(lvard, lvard), obj, method = "BFGS") # numerical estimation
  ssmo <- KFS(mod, smoothing = "state") # state smoothing
  cbind(FLT = ssmo$att[, 1], SMO = ssmo$alphahat[, 1])
}


# Extract the RW component using both the filter and the smoother from a RW + SEAS + WN model.
# It returns a matrix with filtered (FLT) and smoothed (SMO) RW estimates.
trendextractors <- function(x, period = 7) {
  n <- NROW(x)
  lvard <- log(var(x[-1] - x[-n])) # variance of differenced series
  mod <- SSModel(x ~ SSMtrend(1, NA) +
                   SSMseasonal(period, NA), H = NA) # local level model
  obj <- function(pars) { # computes -meanLogLik at parameter values in pars
    mod$Q[1, 1, 1] <<- exp(pars[1])
    mod$Q[2, 2, 1] <<- exp(pars[2])
    mod$H[1, 1, 1] <<- exp(pars[3])
    -logLik(mod) / n
  }
  opt <- optim(c(lvard, -5, lvard), obj, method = "BFGS") # numerical estimation
  ssmo <- KFS(mod, smoothing = "state") # state smoothing
  cbind(FLT = ssmo$att[, 1], SMO = ssmo$alphahat[, 1])
}


# Extract the RW component using both the filter and the smoother from a RW + WN model and omit outliers.
# It returns a matrix with filtered (FLT) and smoothed (SMO) RW estimates.
trendextractor_out <- function(x, th = 2) {
  n <- NROW(x)
  lvard <- log(var(x[-1] - x[-n])) # variance of differenced series
  mod <- SSModel(x ~ SSMtrend(1, NA), H = NA) # local level model
  obj <- function(pars) { # computes -meanLogLik at parameter values in pars
    mod$Q[1, 1, 1] <<- exp(pars[1])
    mod$H[1, 1, 1] <<- exp(pars[2])
    -logLik(mod) / n
  }
  opt <- optim(c(lvard, lvard), obj, method = "BFGS") # numerical estimation
  ssmo <- KFS(mod, smoothing = "disturbance") # state smoothing
  ar <- rstandard(ssmo, "pearson")
  ndx <- which(abs(ar) > th)
  x[ndx] <- NA
  mod <- SSModel(x ~ SSMtrend(1, NA), H = NA) # local level model
  opt <- optim(opt$par, obj, method = "BFGS") # numerical estimation
  ssmo <- KFS(mod, smoothing = "state") # state smoothing
  cbind(FLT = ssmo$att[, 1], SMO = ssmo$alphahat[, 1])
}


# Extract the RW component using both the filter and the smoother from a RW + SEAS + WN model  and omit outliers.
# It returns a matrix with filtered (FLT) and smoothed (SMO) RW estimates.
trendextractors_out <- function(x, period = 7, th = 2) {
  n <- NROW(x)
  lvard <- log(var(x[-1] - x[-n])) # variance of differenced series
  mod <- SSModel(x ~ SSMtrend(1, NA) +
                   SSMseasonal(period, NA), H = NA) # local level model
  obj <- function(pars) { # computes -meanLogLik at parameter values in pars
    mod$Q[1, 1, 1] <<- exp(pars[1])
    mod$Q[2, 2, 1] <<- exp(pars[2])
    mod$H[1, 1, 1] <<- exp(pars[3])
    -logLik(mod) / n
  }
  opt <- optim(c(lvard, -5, lvard), obj, method = "BFGS") # numerical estimation
  ssmo <- KFS(mod, smoothing = c("mean", "disturbance")) # state smoothing
  ar <- rstandard(ssmo, "pearson")
  ndx <- which(abs(ar) > th)
  x[ndx] <- NA
  mod <- SSModel(x ~ SSMtrend(1, NA) +
                   SSMseasonal(period, NA), H = NA) # local level model
  opt <- optim(opt$par, obj, method = "BFGS") # numerical estimation
  ssmo <- KFS(mod, smoothing = "state") # state smoothing
  cbind(FLT = ssmo$att[, 1], SMO = ssmo$alphahat[, 1])
}


# Aggregate time series data by taking fn on m-length chunks
downsample <- function(x, m = 7, fn = mean) {
  n <- NROW(x)
  k <- ceiling(n / m)
  v <- rep(1:k, each = m)[1:n]
  aggregate(x = as.data.frame(x), by = list(v), fn)[, -1]
}


# Extract the RW (trend) component using the smoother from a RW + WN UCM model.
trendextractor_smo <- function(x) {
  n <- NROW(x)
  lvard <- log(var(x[-1] - x[-n])) # variance of differenced series
  mod <- SSModel(x ~ SSMtrend(1, NA), H = NA) # local level model
  obj <- function(pars) { # computes -meanLogLik at parameter values in pars
    mod$Q[1, 1, 1] <<- exp(pars[1])
    mod$H[1, 1, 1] <<- exp(pars[2])
    -logLik(mod) / n
  }
  opt <- optim(c(lvard, lvard), obj, method = "BFGS") # numerical estimation
  ssmo <- KFS(mod, smoothing = "state") # state smoothing
  ssmo$alphahat[, 1]
}


# Extract the RW (trend) component using the filter from a RW + WN UCM model.
trendextractor_flt <- function(x) {
  n <- NROW(x)
  lvard <- log(var(x[-1] - x[-n])) # variance of differenced series
  mod <- SSModel(x ~ SSMtrend(1, NA), H = NA) # local level model
  obj <- function(pars) { # computes -meanLogLik at parameter values in pars
    mod$Q[1, 1, 1] <<- exp(pars[1])
    mod$H[1, 1, 1] <<- exp(pars[2])
    -logLik(mod) / n
  }
  opt <- optim(c(lvard, lvard), obj, method = "BFGS") # numerical estimation
  ssmo <- KFS(mod, smoothing = "none", filtering = "state") # state filter
  ssmo$att[,1]
}


# call trendextractor on every column of x: returns a three dimnesional array
ucm_smo <- function(x) sapply(as.data.frame(x), trendextractor_smo)
ucm_flt <- function(x) sapply(as.data.frame(x), trendextractor_flt)
ucm_two <- function(x) sapply(as.data.frame(x), trendextractor, simplify = "array")

