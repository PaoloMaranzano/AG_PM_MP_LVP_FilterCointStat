library(tidyverse)
library(readr)
library(lubridate)
library(timeDate)
library(readxl)

##### Folder settings
# setwd("C:/Users/paulm/Dropbox/umb visiting programme 2015/Paper 5 - filtering/datasets - vector of 24h per day/Datasets_PUN")

Data <- read_excel("Dataset_electricity_prices_20072021_DE_FR_NL_ES_IT_AT_BE_CH.xlsx")
Data <- Data %>%
  pivot_longer(cols = -dates, names_to = "Series", values_to = "Values") %>% 
  rename(Date = dates) %>%
  mutate(Country = case_when(grepl("de",Series) ~ "DE",
                             grepl("fr",Series) ~ "FR",
                             grepl("it",Series) ~ "IT",
                             grepl("nl",Series) ~ "NL",
                             grepl("ch",Series) ~ "CH",
                             grepl("be",Series) ~ "BE",
                             grepl("es",Series) ~ "ES",
                             grepl("at",Series) ~ "AT"),
         Period = case_when(grepl("deb|frb|itb|nlb|chb|beb|esb|atb",Series) ~ "Base",
                            grepl("dep|frp|itp|nlp|chp|bep|esp|atp",Series) ~ "Peak",
                            grepl("deop1|frop1|itop1|nlop1|chop1|beop1|esop1|atop1",Series) ~ "OffPeak1",
                            grepl("deop2|frop2|itop2|nlop2|chop2|beop2|esop2|atop2",Series) ~ "OffPeak2")) %>%
  select(Date,Country,Period,Values) %>%
  pivot_wider(names_from = Period, values_from = Values) %>%
  rowwise() %>%
  mutate(OffPeak = mean(OffPeak1,OffPeak2,na.rm=T)) %>%
  select(Date, Country, Base, Peak, OffPeak) %>%
  pivot_longer(cols = c(Base,Peak,OffPeak), names_to = "Period", values_to = "Values") %>%
  pivot_wider(names_from = Country, values_from = Values)

