library(KFAS)
library(vars)
library(urca)
library(tidyverse)

setwd("C:/Users/paulm/Dropbox (Personale)/umb visiting programme 2015/Paper 5 - filtering/R")

source("Aux_funs/integration_and_co.R")
source("Aux_funs/functions_for_filtering.R")
source("Aux_funs/functions_for_simulating.R")

results_johansen <- bind_rows(results$johansen)

Sel_rates_k4_r1 <- results_johansen %>%
  filter(k == 4, r == 1) %>%
  group_by(Filter,df,ns_ratio,k,r) %>% 
  group_modify(~ joselect(.,crit = as.numeric(results$crits$`k=4 r=1 df=3 ns_ratio=5`[,3])), .keep = F) %>%
  mutate(df = paste0("df = ",df),
         ns_ratio = paste0("c = ",ns_ratio),
         Filter = factor(x = Filter, levels = c("raw","mean","ucm_flt","ucm_smo"),ordered = T))

Sel_rates_k4_r1 %>%
  dplyr::filter(df %in% c("df = 3", "df = 12"), ns_ratio %in% c("c = 0", "c = 5", "c = 10")) %>%
  ggplot(aes(x = Rels, y = Sel_rate, fill = Filter)) +
  geom_col(position = "dodge") +
  facet_grid(ns_ratio ~ df) +
  labs(y = "Selection rate (%)", x = "Hypothesis") + 
  scale_fill_grey()
