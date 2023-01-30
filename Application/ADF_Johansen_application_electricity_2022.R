library(tidyverse)
library(urca)
library(xts)
library(ggplot2)

server_dismeq <- 0
if (server_dismeq == 1) {
  setwd("C:/Users/paolomaranzano/Desktop/Filtering/Empirical_app_electricity_data/Dataset_2022")
  fld <- "Results2022/"
} else {
  setwd("C:/Users/paulm/Dropbox/umb visiting programme 2015/Paper 5 - filtering/Empirical_app_electricity_data/Dataset_2022") 
  fld <- "Results2022/"
}


source("../../R/Aux_funs/integration_and_co.R")
source("../../R/Aux_funs/functions_for_filtering.R")

estimates <- 1
export_xlsx <- 1
export_output <- 1
outlier <- 1
window_type <- "rolling"
debug <- 0
plots <- 0

periods <- c("Base","Peak","OffPeak")
countries <- c("IT","DE","FR","ES","NL","AT","CH","BE")
countries_ext <- c("Italy","Germany","France","Spain","Netherlands","Austria","Switzerland","Belgium")
if (debug == 1) {
  periods <- c("Base","Peak")
  countries <- c("DE","IT")
  countries_ext <- c("Germany","Italy")
}


ADF_h <- vector("list", length = length(periods))
names(ADF_h) <- c(paste0("ADF_",periods))
Johansen_h <- vector("list", length = length(periods))
names(Johansen_h) <- c(paste0("Johansen_",periods))

if (estimates == 1) {
  for (h in 1:length(periods)) {
    
    source('Data_management_2022.R', encoding = 'UTF-8')
    dt <- Data %>%
      filter(Period == periods[h]) %>%
      select(-Period)
    
    #######################################################
    ########## ADF computation for the countries ##########
    #######################################################
    cat(paste0("ADF for all countries --- Begin at ",Sys.time(),"\n"))
    cnt_idx <- 1:length(countries)
    list_adf <- vector("list", length = length(countries))
    
    for (c in cnt_idx) {
      
      cat(paste0("Country = ",countries[c]," --- Begin at ",Sys.time()))
      
      # lpun <- xts(x = log(dt[,paste0(countries[c],"mean")]), order.by = dt$Date)
      # Gestisco valori negativi mettendoli a mancante (SSM gestiscono missing)
      y <- dt[,countries[c]]
      y[y <= 0] <- NA
      # summary(y)
      lpun <- xts(x = log(y), order.by = dt$Date)
      
      n <- nrow(lpun)
      if (window_type == "rolling") {
        n2 <- 260*3       # Rolling window
      } else {
        n2 <- floor(n/2)  # Recurrent window
      }
      if (debug == 1) {
        n3 <- 10
      } else {
        n3 <- n - n2
      }
      
      ### Stagionalità annuale
      days_year <- 260
      freq365 <- 2 * pi * outer(1:n, 1:16) / days_year
      ### Stagionalità infra-settimanale
      days_week <- 5
      freq5   <- 2 * pi * outer(1:n, 1:2) / days_week
      X <- cbind(cos(freq365), sin(freq365), cos(freq5), sin(freq5))
      
      adfm <- matrix(NA_real_, n3, 4,
                     dimnames = list(time(lpun)[(n2+1):(n2+n3)],
                                     c("RAW", "MEAN", "FLT", "SMO")))
      
      for (t in 1:n3) {
        if (t %% 50 == 0) cat("t =", t, " date =",
                              as.character(time(lpun)[(n2+t)]), "\n")
        reg <- lsfit(X[(1+t):(n2 + t), ], lpun[(1+t):(n2 + t)])
        x <- residuals(reg)
        adfm[t, "RAW"]  <- adf_test(x, "const", days_week*3, "AICC")$adf
        adfm[t, "MEAN"] <- adf_test(downsample(x), "const", 1*3, "AICC")$adf
        if (outlier == 1) {
          fs <- trendextractors_out(x)
        } else {
          fs <- trendextractors(x) 
        }
        adfm[t, "FLT"]  <- adf_test(fs[, "FLT"], "const", days_week*3, "AICC")$adf
        adfm[t, "SMO"]  <- adf_test(fs[, "SMO"], "const", days_week*3, "AICC")$adf
      }
      
      adf_tidy <- adfm %>%
        data.frame() %>%
        mutate(Date = dt$Date[(n2+1):(n2+n3)]) %>%
        pivot_longer(cols = 1:4, names_to = "Filter", values_to = "ADF_Tau_stat") %>%
        mutate(Country = countries_ext[c])
      
      list_adf[[c]] <- adf_tidy
      
      if (export_xlsx == 1) {
        paste0(fld,paste0(window_type,"_ADF_",countries[c],"_",Sys.Date(),".csv"))
      }
      
      cat(paste0("Country = ",countries[c]," --- Ended at ",Sys.time()))
    }
    ### Bind ADF results for all the countries
    adf_tidy_full <- bind_rows(list_adf)
    adf_tidy_full <- adf_tidy_full %>%
      mutate(Filter = recode(Filter,"RAW" = "no","SMO" = "ucm_smo","FLT"="ucm_flt","MEAN"="mean"))
    ### Store ADF results for h-th hour
    ADF_h[[h]] <- adf_tidy_full
    cat(paste0("ADF for all countries --- Ended at ",Sys.time(),"\n"))
    
    
    
    ############################################################
    ########## Johansen computation for the countries ##########
    ############################################################
    cat(paste0("Johansen for all countries --- Begin at ",Sys.time(),"\n"))
    
    ym <- dt[,countries]
    ym[ym <= 0] <- NA
    lpun <- xts(x = log(ym), order.by = dt$Date)
    
    n <- nrow(lpun)
    if (window_type == "rolling") {
      n2 <- 260*3       # Rolling window
    } else {
      n2 <- floor(n/2)  # Recurrent window
    }
    if (debug == 1) {
      n3 <- 2
    } else {
      n3 <- n - n2
    }
    
    freq365 <- 2 * pi * outer(1:n, 1:16) / 260
    freq5   <- 2 * pi * outer(1:n, 1:2) / days_week
    X <- cbind(cos(freq365), sin(freq365), cos(freq5), sin(freq5))
    
    max_num_coint_rels <- length(countries)-1
    
    joh_raw <- matrix(NA_real_, n3, max_num_coint_rels+1, dimnames = list(time(lpun)[(n2+1):(n2+n3)],paste0("k",0:max_num_coint_rels)))
    joh_mean <- matrix(NA_real_, n3, max_num_coint_rels+1, dimnames = list(time(lpun)[(n2+1):(n2+n3)],paste0("k",0:max_num_coint_rels)))
    joh_flt <- matrix(NA_real_, n3, max_num_coint_rels+1, dimnames = list(time(lpun)[(n2+1):(n2+n3)],paste0("k",0:max_num_coint_rels)))
    joh_smo <- matrix(NA_real_, n3, max_num_coint_rels+1, dimnames = list(time(lpun)[(n2+1):(n2+n3)],paste0("k",0:max_num_coint_rels)))
    
    for (t in 1:n3) {
      if (t %% 50 == 0) cat("iteration =", t, " date =", as.character(time(lpun)[(n2+t)]), "\n")
      reg <- lsfit(X[(1+t):(n2 + t), ], lpun[(1+t):(n2 + t)])
      x <- residuals(reg)
      x_raw <- x
      x_mean <- downsample(x)
      if (outlier == 1) {
        x_flt <- apply(x, 2, function(x) trendextractors_out(x)[,"FLT"])
        x_smo <- apply(x, 2, function(x) trendextractors_out(x)[,"SMO"])
      } else {
        x_flt <- apply(x, 2, function(x) trendextractors(x)[,"FLT"])
        x_smo <- apply(x, 2, function(x) trendextractors(x)[,"SMO"]) 
      }
      joh_raw[t,]  <- johansen_test(y = x, type = "uconst", lags = 5*3)$statistics[,"trace"]
      joh_mean[t,]  <- johansen_test(y = x_mean, type = "uconst", lags = 3)$statistics[,"trace"]
      joh_flt[t,]  <- johansen_test(y = x_flt, type = "uconst", lags = 5*3)$statistics[,"trace"]
      joh_smo[t,]  <- johansen_test(y = x_smo, type = "uconst", lags = 5*3)$statistics[,"trace"]
    }
    
    joh_raw <- joh_raw %>%
      data.frame() %>%
      mutate(Date = dt$Date[(n2+1):(n2+n3)], Filter = "Raw")
    joh_mean <- joh_mean %>%
      data.frame() %>%
      mutate(Date = dt$Date[(n2+1):(n2+n3)], Filter = "Mean")
    joh_flt <- joh_flt %>%
      data.frame() %>%
      mutate(Date = dt$Date[(n2+1):(n2+n3)], Filter = "KF")
    joh_smo <- joh_smo %>%
      data.frame() %>%
      mutate(Date = dt$Date[(n2+1):(n2+n3)], Filter = "KS")
    ### Bind Johansen results for all the countries
    johansen <- bind_rows(joh_raw,joh_mean,joh_flt,joh_smo)
    johansen_tidy <- johansen %>%
      mutate(Filter = recode(Filter,"Raw" = "no","KS" = "ucm_smo","KF"="ucm_flt","Mean"="mean")) %>%
      pivot_longer(cols = 1:length(countries), names_to = "K") %>%
      mutate(K = recode(K,"k0"="r = 0","k1"="r = 1","k2"="r = 2","k3"="r = 3","k4"="r = 4"))
    ### Store Johansen results for h-th hour
    Johansen_h[[h]] <- johansen_tidy
    cat(paste0("Johansen for all countries --- Ended at ",Sys.time(),"\n"))
  }
} else {
  load(paste0(fld,"c.RData"))
}


#################################
########## Save output ##########
#################################

if (export_output == 1) {
  save.image(paste0(fld,paste0("Res_",window_type,"_",paste0(countries,collapse = "_"),"_",Sys.Date(),".RData")))
}



