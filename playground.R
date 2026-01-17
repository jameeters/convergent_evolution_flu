library(ape)
library(dplyr)
library(lubridate)
library(stringr)

setwd("/Users/james/Documents/convergent/out")

treetime = read.nexus('treetime/timetree.nexus')


dates = read.csv('treetime/dates.tsv', sep='\t')
dates = dates %>% 
  rename(node = X.node) %>% 
  mutate(date = as.Date(ymd(date)))

plot(treetime, show.tip.label=F)

