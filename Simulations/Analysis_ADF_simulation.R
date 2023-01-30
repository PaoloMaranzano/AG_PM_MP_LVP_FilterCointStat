library(KFAS)
library(vars)
library(urca)
library(tidyverse)
library(ggplot2)

setwd("C:/Users/paulm/Dropbox (Personale)/umb visiting programme 2015/Paper 5 - filtering/R")

source("Aux_funs/integration_and_co.R")
source("Aux_funs/functions_for_filtering.R")
source("Aux_funs/functions_for_simulating.R")

########## Size
results_ADF <- read.csv("sim_adf_tests_H0.csv")
results_ADF$df <- paste0("df=", results_ADF$df)
results_ADF$df <- ordered(x = results_ADF$df, levels = c("df=3","df=6","df=9","df=12"))

pdf(file = paste0("graph_adf_rejections.pdf"), width = 8, height = 6)
ggplot(results_ADF, aes(x = ns_ratio, y = rejection_rate, color = filter, linetype = filter)) +
  geom_line(size = 1) +
  facet_wrap(~df, 2, 2) +
  xlab("noise to signal ratio") +
  ylab("rejection rate") +
  scale_x_continuous(breaks = 0:10, minor_breaks = NULL) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), minor_breaks = seq(0, 1, 0.05))
dev.off()

########## Size-adjusted power
results_ADF <- read.csv("sim_adf_tests_H1.csv")
results_ADF$df <- paste0("df=", results_ADF$df)
results_ADF$df <- ordered(x = results_ADF$df, levels = c("df=3","df=6","df=9","df=12"))

pdf(file = paste0("graph_adf_powers.pdf"), width = 8, height = 6)
ggplot(results_ADF, aes(x = ns_ratio, y = rejection_rate, color = filter, linetype = filter)) +
  geom_line(size = 1) +
  facet_wrap(~df, 2, 2) +
  xlab("noise to signal ratio") +
  ylab("rejection rate") +
  scale_x_continuous(breaks = 0:10, minor_breaks = NULL) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), minor_breaks = seq(0, 1, 0.1))
dev.off()
