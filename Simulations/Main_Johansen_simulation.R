library(KFAS)
library(vars)
library(urca)
library(tidyverse)

setwd("C:/Users/paulm/Dropbox (Personale)/umb visiting programme 2015/Paper 5 - filtering/R/Simulations/")

source("../Aux_funs/integration_and_co.R")
source("../Aux_funs/functions_for_filtering.R")
source("../Aux_funs/functions_for_simulating.R")


##### Simulation scheme 1: one cointegrating relationship
beta  <- cbind(b1 = c(1, -1, 1, -1))
alpha <- cbind(a1 = c(0, 0.1, -0.1, 0.1))
Pi <- alpha %*% t(beta)
phi <- diag(4) + Pi
eigen(phi)


##### Test per H relazioni di cointegrazione #####
## Fisso il numero di simulazioni da effettuare (sim)
nsim <- 10000
n <- 365 * 3
k <- 4
r <- 1
vdf <- seq(3, 12, 3)
# vdf <- 3
vc <- 0:10
# vc <- 0

##### Simulations
results <- list()
set.seed(123)
for (df in vdf) {
  for (ns_ratio in vc) {
    cat("df =", df, " ns_ratio =", ns_ratio, "--", as.character(Sys.time()), "\n")
    X <- sim_noisyVAR(phi = phi, nsim, n, ns_ratio, df, k, r)
    
    joh_raw <- joh_mean <- joh_flt <- joh_smo <- matrix(NA_real_, nsim, k, 
                                                        dimnames = list(1:nsim,paste0("k_star",0:(k-1))))
    
    # Critical from A1 by Johansen&Juselius (1990) and Osterwald-Lenum (1992, pag. 8)
    # unconstrained constant
    Crit_Joh <- data.frame(c10 = c(150.53,118.50,89.48,64.84,43.95,26.79,13.33,2.69),
                           c5 = c(156.00,124.24,94.15,68.52,47.21,29.69,15.41,3.76),
                           c1 = c(168.36,133.57,103.18,76.07,54.46,35.65,20.04,6.65))
    Crit_Joh <- Crit_Joh[seq(dim(Crit_Joh)[1]-(k-1),dim(Crit_Joh)[1]),]
    K <- paste("r =",0:(k-1))
    crits <- bind_cols(data.frame(K),Crit_Joh)
    
    for (s in 1:nsim) {
      # Filtering series
      X_raw <- X[[s]]
      X_mean <- downsample(X[[s]])
      X_flt <- apply(X[[s]], 2, function(x) trendextractor(x)[,"FLT"])
      X_smo <- apply(X[[s]], 2, function(x) trendextractor(x)[,"SMO"])
      # Store results from Johansen
      joh_raw[s,]  <- johansen_test(y = X_raw, type = "uconst", lags = vars::VARselect(X_raw, 20, "none")$selection[1])$statistics[,"trace"]
      joh_mean[s,]  <- johansen_test(y = X_mean, type = "uconst", lags = vars::VARselect(X_mean, 20, "none")$selection[1])$statistics[,"trace"]
      joh_flt[s,]  <- johansen_test(y = X_flt, type = "uconst", lags = vars::VARselect(X_flt, 20, "none")$selection[1])$statistics[,"trace"]
      joh_smo[s,]  <- johansen_test(y = X_smo, type = "uconst", lags = vars::VARselect(X_smo, 20, "none")$selection[1])$statistics[,"trace"]
    }
    
    joh_raw <- joh_raw %>%
      data.frame() %>%
      mutate(Filter = "raw",df = df, ns_ratio = ns_ratio, k = k, r = r)
    joh_mean <- joh_mean %>%
      data.frame() %>%
      mutate(Filter = "mean",df = df, ns_ratio = ns_ratio, k = k, r = r)
    joh_flt <- joh_flt %>%
      data.frame() %>%
      mutate(Filter = "ucm_flt",df = df, ns_ratio = ns_ratio, k = k, r = r)
    joh_smo <- joh_smo %>%
      data.frame() %>%
      mutate(Filter = "ucm_smo",df = df, ns_ratio = ns_ratio, k = k, r = r)
    
    johansen <- bind_rows(joh_raw,joh_mean,joh_flt,joh_smo)
    rownames(johansen) <- NULL
    
    results$johansen[[paste0("k=",k," r=",r," df=", df, " ns_ratio=", ns_ratio)]] <- johansen
    results$crits[[paste0("k=",k," r=",r," df=", df, " ns_ratio=", ns_ratio)]] <- crits
  }
}

##### Save output
save(results,beta,alpha,phi,nsim,n,file = paste0("Sim_Johansen_k",k,"_r",r,".RData"))