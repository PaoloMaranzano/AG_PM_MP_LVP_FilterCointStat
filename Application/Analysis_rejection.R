library(tidyverse)
library(ggplot2)
library(xtable)

setwd("C:/Users/paulm/Dropbox/umb visiting programme 2015/Paper 5 - filtering/Empirical_app_electricity_data/Dataset_2022")
fld <- "Results2022/"

# load("Results2022/Res_rolling_IT_DE_FR_ES_NL_AT_CH_BE_2022-03-04.RData")
load("Results2022/Res_rolling_DE_FR_NL_AT_2022-03-10.RData")

plots <- 0


############################################
########## ADF and Johansen plots ##########
############################################
if (plots == 1) {
  ##### Colors
  cols_flt <- c("ucm_smo" = "#0000FF", "ucm_flt" = "#00FF33", "mean" = "#FF9933", "no" = "#666666")
  
  for (h in 1:length(periods)) {
    ##### Plot ADF
    period <- periods[h]
    plot_adf <- ADF_h[[paste0("ADF_",period)]] %>%
      data.frame() %>%
      # filter(Country == "Germany") %>%
      ggplot(aes(x = Date, y = ADF_Tau_stat, color = Filter)) +
      geom_line() + 
      facet_wrap(~ Country,scales = "free") + 
      labs(x = "Date", y = "ADF Tau statistics",
           title = paste0("Rolling window ADF statistics (",period,")"),
           subtitle = paste(countries_ext,collapse = " - ")) +
      theme(legend.position = "bottom",plot.title = element_text(size = 18)) + 
      scale_color_manual(values = cols_flt) + 
      geom_hline(yintercept = c(-3.435466,    # Critico ADF al 1%
                                -2.862951,    # Critico ADF al 5%
                                -2.567547),   # Critico ADF al 10%
                 lwd=1.1)
    pdf(file = paste0(fld,"ADF_",paste(countries,collapse = "_"),"_",period,".pdf"),width = 10, height = 8)
    print(plot_adf)
    dev.off()
    
    
    ##### Plot Johansen
    # Critical from A1 by Johansen&Juselius (1990) and Osterwald-Lenum (1992, pag. 8)
    Crit_Joh <- data.frame(c10 = c(150.53,118.50,89.48,64.84,43.95,26.79,13.33,2.69),
                           c5 = c(156.00,124.24,94.15,68.52,47.21,29.69,15.41,3.76),
                           c1 = c(168.36,133.57,103.18,76.07,54.46,35.65,20.04,6.65))
    Crit_Joh <- Crit_Joh[seq(dim(Crit_Joh)[1]-(length(countries)-1),dim(Crit_Joh)[1]),]
    K <- paste("r =",0:(length(countries)-1))
    Crit_Joh <- bind_cols(data.frame(K),Crit_Joh)
    
    plot_johansen <- Johansen_h[[paste0("Johansen_",period)]] %>%
      mutate(K = case_when(K == "k5" ~ "r = 5",
                           K == "k6" ~ "r = 6",
                           K == "k7" ~ "r = 7",
                           TRUE ~ K)) %>%
      ggplot(aes(x = Date, y = value, color = Filter)) +
      geom_line() + 
      labs(x = "Date", y = "Trace statistics",
           title = paste0("Rolling window Johansen trace statistics (",period,")"),
           subtitle = paste(countries_ext,collapse = " - ")) +
      theme(legend.position = "bottom",plot.title = element_text(size = 18)) + 
      scale_color_manual(values = cols_flt) + 
      geom_hline(aes(yintercept = c10), Crit_Joh, lwd=1.1) + 
      geom_hline(aes(yintercept = c5), Crit_Joh, lwd=1.1) + 
      geom_hline(aes(yintercept = c1), Crit_Joh, lwd=1.1) + 
      facet_wrap(~ K, scales = "free")
    pdf(file = paste0(fld,"Johansen_",paste(countries,collapse = "_"),"_",period,".pdf"), width = 10, height = 8)
    print(plot_johansen)
    dev.off()
  }
  
  
  #####################################################
  ########## Grafico serie storica originale ##########
  #####################################################
  plot_dt <- dt %>%
    pivot_longer(cols = -1, names_to = "Country", values_to = "Value") %>%
    filter(Country %in% countries) %>%
    mutate(Country = recode(Country,
                            "AT" = "Austria",
                            "DE"="Germany",
                            "ES"="Spain",
                            "FR"="France",
                            "IT"="Italy",
                            "CH" = "Switzerland",
                            "BE" = "Belgium",
                            "NL"="Netherland")) %>%
    ggplot(aes(x = Date, y = Value)) +
    geom_line() + 
    facet_wrap(~ Country,nrow = 3, scales = "free") + 
    labs(x = "", y = "???/MWh", title = "") + 
    theme_minimal() + 
    theme(strip.text = element_text(size=14),axis.title.y = element_text(size=12))
  pdf(file = paste0(fld,"TS_prices_",paste(countries,collapse = "_"),".pdf"),width = 10, height = 8)
  print(plot_dt)
  dev.off()
}







h <- 1

for (h in 1:length(periods)) {
  ##### ADF test for each scenario
  period <- periods[h]
  adf_tidy_full <- ADF_h[[paste0("ADF_",period)]] %>%
    mutate(Filter = recode(Filter,
                           "no"="No filters",
                           "mean"="Moving average",
                           "ucm_flt"="Kalman filter",
                           "ucm_smo"="Kalman smoother"))

  ###################################################
  ########## ADF tests for the full sample ##########
  ###################################################
  adf_tidy_full <- adf_tidy_full %>%
    #  filter(Date >= "2011-03-01", Date <= "2014-09-01") %>%
    mutate(Reject10 = ifelse(ADF_Tau_stat < -2.567547,1,0),
           Reject5 = ifelse(ADF_Tau_stat < -2.862951,1,0),
           Reject1 = ifelse(ADF_Tau_stat < -3.435466,1,0))
  ##### Rejection share (5%) by country and filter
  ADF_5 <- adf_tidy_full %>%
    group_by(Country,Filter) %>%
    summarise(n = n(),
              Rej_share5 = sum(Reject5)/n*100) %>%
    pivot_wider(names_from = Country, values_from = Rej_share5)
  ADF_fullsample_5_latex <- xtable(ADF_5,
                                   caption = paste0("ADF test on the full sample (1st January 2010 - 31th December 2021): rejection rates at 5% (",paste(countries,collapse = "-"),") - ",period), 
                                   align=c("l","l","l",rep("c",length(countries))))
  ADF_fullsample_5_latex
  
  ########################################################
  ########## ADF tests for the 2011-2014 period ##########
  ########################################################
  adf_tidy_red <- adf_tidy_full %>%
    filter(Date >= "2011-03-01", Date <= "2014-09-01") %>%
    mutate(Reject10 = ifelse(ADF_Tau_stat < -2.567547,1,0),
           Reject5 = ifelse(ADF_Tau_stat < -2.862951,1,0),
           Reject1 = ifelse(ADF_Tau_stat < -3.435466,1,0))
  ##### Rejection share (5%) by country and filter
  ADF_5 <- adf_tidy_red %>%
    group_by(Country,Filter) %>%
    summarise(n = n(),
              Rej_share5 = sum(Reject5)/n*100) %>%
    pivot_wider(names_from = Country, values_from = Rej_share5)
  ADF_redsample_5_latex <- xtable(ADF_5,
                                  caption = paste0("ADF test on the full sample (1st January 2010 - 31th December 2021): rejection rates at 5% (",paste(countries,collapse = "-"),") - ",period), 
                                  align=c("l","l","l",rep("c",length(countries))))
  ADF_redsample_5_latex
  
  
  
  ######################################
  ########## Johansen's tests ##########
  ######################################
  # Critical from A1 by Johansen&Juselius (1990) and Osterwald-Lenum (1992, pag. 8)
  Crit_Joh <- data.frame(c10 = c(150.53,118.50,89.48,64.84,43.95,26.79,13.33,2.69),
                         c5 = c(156.00,124.24,94.15,68.52,47.21,29.69,15.41,3.76),
                         c1 = c(168.36,133.57,103.18,76.07,54.46,35.65,20.04,6.65))
  Crit_Joh <- Crit_Joh[seq(dim(Crit_Joh)[1]-(length(countries)-1),dim(Crit_Joh)[1]),]
  K <- paste("r =",0:(length(countries)-1))
  Crit_Joh <- bind_cols(data.frame(K),Crit_Joh)
  
  ##### Johansen test for each scenario
  johansen_tidy_full <- Johansen_h[[paste0("Johansen_",period)]] %>%
    mutate(K = recode(K,
                      "r=0" = "r = 0",
                      "r=1" = "r = 1",
                      "r=2" = "r = 2",
                      "r=3" = "r = 3",
                      "r=4" = "r = 4"),
           Filter = recode(Filter,
                           "no"="No filters",
                           "mean"="Moving average",
                           "ucm_flt"="Kalman filter",
                           "ucm_smo"="Kalman smoother"))
  
  johansen_tidy_full <- Johansen_h[[paste0("Johansen_",period)]] %>%
    mutate(K = case_when(K == "k5" ~ "r = 5",
                         K == "k6" ~ "r = 6",
                         K == "k7" ~ "r = 7",
                         TRUE ~ K),
           K = recode(K,
                      "r=0" = "r = 0",
                      "r=1" = "r = 1",
                      "r=2" = "r = 2",
                      "r=3" = "r = 3",
                      "r=4" = "r = 4"),
           Filter = recode(Filter,
                           "no"="No filters",
                           "mean"="Moving average",
                           "ucm_flt"="Kalman filter",
                           "ucm_smo"="Kalman smoother"))
  
  johansen_tidy_full <- left_join(johansen_tidy_full,Crit_Joh, by = "K")
  
  johansen_tidy_full <- johansen_tidy_full %>%
    mutate(Reject10 = ifelse(value > c10,1,0),
           Reject5 = ifelse(value > c5,1,0),
           Reject1 = ifelse(value > c1,1,0))
  
  col_filters <- c("No filters" = "#FF6600", "Moving average" = "#33CC00",
                   "Kalman filter" = "#FF0000", "Kalman smoother" = "#0000FF")
  

  p_joh <- johansen_tidy_full %>%
    mutate(Order = substr(K,5,5)) %>%
    select(Date,Filter,Reject1,Order) %>%
    filter(Reject1 == 0) %>%
    group_by(Date,Filter) %>%
    summarise(min_ord = as.numeric(min(Order))) %>%
    ggplot(mapping = aes(x = Date, y = min_ord, col=Filter)) +
    geom_point() + 
    geom_smooth(col="black",size=2) + 
    facet_wrap( ~ Filter,nrow = 4) + 
    scale_color_manual(values = col_filters) + 
    theme(legend.position = "") + 
    theme(axis.text.x=element_text(angle=45, hjust=1)) + 
    scale_x_datetime(breaks= lubridate::ymd_hms(paste0(c(2010:2022),"-01-01 00:00:00"))) + 
    labs(title = paste0(period," hour: ",paste(countries_ext,collapse = "-")),y="Cointegration relationships",x = "")
    # labs(title = paste0("Number of cointegration relationships identified by the filters (",paste(countries,collapse = "-"),")"),
    #      y="Cointegration relationships",x = "")
  pdf(file = paste0("Number_coint_rels_identified_",paste(countries,collapse = "_"),"_",period,".pdf"), width = 10, height = 8)
  print(p_joh)
  dev.off()
  
}


