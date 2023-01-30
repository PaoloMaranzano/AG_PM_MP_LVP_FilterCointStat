library(KFAS)
library(vars)
library(urca)
library(tidyverse)

setwd("C:/Users/paulm/Dropbox (Personale)/umb visiting programme 2015/Paper 5 - filtering/R/Simulations/")

source("../Aux_funs/integration_and_co.R")
source("../Aux_funs/functions_for_filtering.R")
source("../Aux_funs/functions_for_simulating.R")


####################################################
########## ADF test for noisy time series ##########
####################################################
nsim <- 100
n <- 365 * 3
# vdf <- seq(3, 12, 3)
vdf <- 3
vc <- 0:10

##### Simulations under H0 (actual size) #####
results <- list()
set.seed(123)
for (df in vdf) {
  for (ns_ratio in vc) {
    cat("df =", df, " ns_ratio =", ns_ratio, "--", as.character(Sys.time()), "\n")
    X <- sim_noisyrw(nsim = nsim, n = n, df = df, ns_ratio = ns_ratio)
    
    ADF_stat <- matrix(NA_real_, nsim, 4, dimnames = list(1:nsim,c("raw","mean","ucm_flt","ucm_smo")))
    
    # unconstrained constant
    crits <- dfcv(n = n, type = "none")
    names(crits) <- c("c1","c5","c10")
    
    # Filtering series
    X_raw <- X
    X_mean <- downsample(X)
    X_flt <- apply(X, 2, function(x) trendextractor(x)[,"FLT"])
    X_smo <- apply(X, 2, function(x) trendextractor(x)[,"SMO"])
    # Store results
    ADF_stat[,1]  <- apply(X_raw, 2, function(x) adf_test(x, type = "none", lags = 10, selectlags = "AICC")$adf)
    ADF_stat[,2]  <- apply(X_mean, 2, function(x) adf_test(x, type = "none", lags = 10, selectlags = "AICC")$adf)
    ADF_stat[,3]  <- apply(X_flt, 2, function(x) adf_test(x, type = "none", lags = 10, selectlags = "AICC")$adf)
    ADF_stat[,4]  <- apply(X_smo, 2, function(x) adf_test(x, type = "none", lags = 10, selectlags = "AICC")$adf)
    
    ADF_stat <- ADF_stat %>%
      data.frame() %>%
      mutate(df = df, ns_ratio = ns_ratio)
    
    results$ADF[[paste0("df=", df, " ns_ratio=", ns_ratio)]] <- ADF_stat
    results$crits[[paste0("df=", df, " ns_ratio=", ns_ratio)]] <- crits
  }
}
# Save output
save(results,nsim,n,file = paste0("Sim_ADF_H0.RData"))
# Aggregation to rejection rates
results_ADF <- bind_rows(results$ADF)
Rej_rates_ADF <- results_ADF %>%
  pivot_longer(cols = c(raw,mean,ucm_flt,ucm_smo), names_to = "Filter", values_to = "Tau") %>%
  group_by(Filter,df,ns_ratio) %>% 
  summarise(Rej_rate = mean(Tau < results$crits[[1]][2])*100) %>%
  mutate(df = paste0("df = ",df),
         Filter = factor(x = Filter, levels = c("raw","mean","ucm_flt","ucm_smo"),ordered = T))



##### Simulations under H1 (size-adjusted power) #####
results <- list()
set.seed(123)
for (df in vdf) {
  for (ns_ratio in vc) {
    cat("df =", df, " ns_ratio =", ns_ratio, "--", as.character(Sys.time()), "\n")
    X <- sim_noisyAR1(nsim = nsim, n = n, df = df, ns_ratio = ns_ratio, phi = 0.98)
    
    ADF_stat <- matrix(NA_real_, nsim, 4, dimnames = list(1:nsim,c("raw","mean","ucm_flt","ucm_smo")))
    
    # unconstrained constant
    crits <- dfcv(n = n, type = "none")
    names(crits) <- c("c1","c5","c10")
    
    # Filtering series
    X_raw <- X
    X_mean <- downsample(X)
    X_flt <- apply(X, 2, function(x) trendextractor(x)[,"FLT"])
    X_smo <- apply(X, 2, function(x) trendextractor(x)[,"SMO"])
    # Store results
    ADF_stat[,1]  <- apply(X_raw, 2, function(x) adf_test(x, type = "none", selectlags = "AICC")$adf)
    ADF_stat[,2]  <- apply(X_mean, 2, function(x) adf_test(x, type = "none", selectlags = "AICC")$adf)
    ADF_stat[,3]  <- apply(X_flt, 2, function(x) adf_test(x, type = "none", selectlags = "AICC")$adf)
    ADF_stat[,4]  <- apply(X_smo, 2, function(x) adf_test(x, type = "none", selectlags = "AICC")$adf)
    
    ADF_stat <- ADF_stat %>%
      data.frame() %>%
      mutate(df = df, ns_ratio = ns_ratio)
    
    results$ADF[[paste0("df=", df, " ns_ratio=", ns_ratio)]] <- ADF_stat
    results$crits[[paste0("df=", df, " ns_ratio=", ns_ratio)]] <- crits
  }
}
# Save output
save(results,nsim,n,file = paste0("Sim_ADF_H1.RData"))
# Aggregation to rejection rates
results_ADF <- bind_rows(results$ADF)
Rej_rates_ADF <- results_ADF %>%
  pivot_longer(cols = c(raw,mean,ucm_flt,ucm_smo), names_to = "Filter", values_to = "Tau") %>%
  group_by(Filter,df,ns_ratio) %>% 
  summarise(Rej_rate = mean(Tau < results$crits[[1]][2])*100) %>%
  mutate(df = paste0("df = ",df),
         Filter = factor(x = Filter, levels = c("raw","mean","ucm_flt","ucm_smo"),ordered = T))

