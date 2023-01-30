# simulate nsim time series of length n of type RW + ns_ratio * student(df)
sim_noisyrw <- function(nsim, n, df, ns_ratio) {
  eps <- matrix(rnorm(nsim * n), n, nsim)
  timeSeries::colCumsums(eps) + ns_ratio * matrix(rt(nsim * n, df = df), n, nsim) * sqrt((df - 2) / df)
}


# simulate nsim time series of length n of type AR(1) + ns_ratio * student(df)
sim_noisyAR1 <- function(nsim, n, df, ns_ratio, phi) {
  eps <- matrix(rnorm(nsim * n), n, nsim)
  stats::filter(x = eps, filter = phi, method = "recursive") + ns_ratio * matrix(rt(nsim * n, df = df), n, nsim) * sqrt((df - 2) / df)
}


### Function for simulating nsim noisy VAR + ns_ratio * student(df)
#   phi: either matrix of VAR(1) or NULL to simulate it
#   nsim: number of simulations
#   n: length of time series
#   c: noise / signal variance ratio
#   df: degrees of freedom of Student's t
#   flt: function to filter time series
#   k: number of time series
#   r: number of coitegrating relations
#   df = gdl disturbo leptocurtico = 3, 6, 12
#   c = 1, 2, 4, 8, 16
# genero rumore t student e lo standardizzo (divido per sigma(t))
sim_noisyVAR <- function(phi = NULL, nsim, n, df, ns_ratio, k, r) {
  flag  <- is.null(phi)
  if (!flag) k <- dim(phi)[1]
  
  Y <- vector(mode = "list", length = nsim)
  for (s in 1:nsim) {
    Y[[s]] <- matrix(0, nrow = n, ncol = k)
  }
  
  # X <- matrix(rnorm(n * k), n, k, dimnames = list(NULL, paste0("Y", 1:k)))
  tsd <- sqrt(df / (df - 2)) # st.dev of Student's t

  for (s in 1:nsim) {
    if (flag) { # radom VAR(1) matrix with 1's and 0.5's eigenvalues
      PI     <- matrix(runif(k * r), k, r) %*% matrix(runif(k * r), r, k)
      phi <- diag(k) + PI
      eigenv <- c(rep(1, k - r), rep(0.5, r))
      phi <- Re(eigen(phi)$vectors %*% diag(eigenv) %*% solve(eigen(phi)$vectors))
    }
    # Y1 <- VAR.sim(B = phi, n = n, include = "none")
    Y1 <- varFilter(t(phi), matrix(rnorm(n * k), n, k), matrix(0, 1, k))[-1,]
    alpha <- ns_ratio * apply(diff(as.matrix(Y1)), MARGIN=2, FUN=var)
    noise <- matrix(rt(df = df, n = n * k), nrow = n, ncol = k)
    # Y1 <- Y1 + rep(alpha, each = n) * noise / tsd
    dimnames(Y1) <- list(NULL, paste0("Y", 1:k))
    Y[[s]] <- Y1 + rep(alpha, each = n) * noise / tsd
  }

  return(Y)
}




### Function for simulating Johansen test
#   phi: either matrix of VAR(1) or NULL to simulate it
#   nsim: number of simulations
#   n: length of time series
#   c: noise / signal variance ratio
#   df: degrees of freedom of Student's t
#   flt: function to filter time series
#   k: number of time series
#   r: number of coitegrating relations
#   df = gdl disturbo leptocurtico = 3, 6, 12
#   c = 1, 2, 4, 8, 16
# genero rumore t student e lo standardizzo (divido per sigma(t))
test <- function(phi = NULL, nsim, n, ns_ratio, df, flt1, flt2, k, r, ...) {
  flag  <- is.null(phi)
  if (!flag) k <- dim(phi)[1]
  
  stats <- vector(mode = "list", length = length(flts))
  for (nf in 1:length(stats)) {
    stats[[nf]] <- matrix(0, nrow = nsim, ncol = k)
  }
  
  # stats <- matrix(0, nrow = nsim, ncol = k)
  # stats.filt1 <- matrix(0, nrow = nsim, ncol = k)
  # stats.filt2 <- matrix(0, nrow = nsim, ncol = k)
  
  X <- matrix(rnorm(n * k), n, k, dimnames = list(NULL, paste0("Y", 1:k)))
  c1  <- as.numeric(ca.jo(X, type="trace")@cval[, 3])
  c5  <- as.numeric(ca.jo(X, type="trace")@cval[, 2])
  c10 <- as.numeric(ca.jo(X, type="trace")@cval[, 1])
  
  tsd <- sqrt(df / (df - 2)) # st.dev of Student's t
  for (row in 1:nsim) {
    if (flag) { # radom VAR(1) matrix with 1's and 0.5's eigenvalues
      PI     <- matrix(runif(k * r), k, r) %*% matrix(runif(k * r), r, k)
      phi <- diag(k) + PI
      eigenv <- c(rep(1, k - r), rep(0.5, r))
      phi <- Re(eigen(phi)$vectors %*% diag(eigenv) %*% solve(eigen(phi)$vectors))
    }
    # Y1 <- VAR.sim(B = phi, n = n, include = "none")
    Y1 <- varFilter(t(phi), matrix(rnorm(n * k), n, k), matrix(0, 1, k))[-1,]
    alpha <- ns_ratio * apply(diff(as.matrix(Y1)), MARGIN=2, FUN=var)
    noise <- matrix(rt(df = df, n = n * k), nrow = n, ncol = k)
    # Y <- as.data.frame(Y1 + rep(alpha, each = n) * noise / tsd)
    Y <- Y1 + rep(alpha, each = n) * noise / tsd
    # stats[row, ] <- ca.jo(x = Y, type = "trace", ecdet = "none")@teststat
    lags <- VARselect(Y, 10, "none")$selection[1]
    
    # for (nf in 1:length(stats)) {
    #   stats[[nf]][row,] <- johansen_test(Y, "uconst", lags)$statistics[, 2]
    # }
    # stats[row, ] <- johansen_test(Y, "uconst", lags)$statistics[, 2]
    # Y.filt <- apply(Y, 2,  FUN = flt, ...)
    # Y.filt1 <- flt1(Y, ...)
    # stats.filt1[row, ] <- johansen_test(Y.filt1, "uconst", lags)$statistics[, 2]
    # Y.filt2 <- flt2(Y, ...)
    # stats.filt2[row, ] <- johansen_test(Y.filt2, "uconst", lags)$statistics[, 2]
  }
  stats <- data.frame(stats)
  # colnames(stats) <- paste0("r", (k - 1):0)
  colnames(stats) <- paste0("r", 0:(k - 1))
  
  # stats.filt <- data.frame(stats.filt)
  # colnames(stats.filt) <- paste0("r", (k - 1):0)
  stats.filt1 <- data.frame(stats.filt1)
  colnames(stats.filt1) <- paste0("r", 0:(k - 1))
  stats.filt2 <- data.frame(stats.filt2)
  colnames(stats.filt2) <- paste0("r", 0:(k - 1))
  
  crit <- data.frame(rbind(t(c1[k:1]), t(c5[k:1]), t(c10[k:1])))
  colnames(crit) <- paste0("c.r", 0:(k - 1))
  
  return(list(stats = stats,
              stats.filt1 = stats.filt1,
              stats.filt2 = stats.filt2,
              crit=crit))
}