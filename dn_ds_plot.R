library(dplyr)
library(ggplot2)
library(ggrepel)

theme_set(theme_bw())
setwd('/Users/james/Documents/convergent/out')

mut_counts = read.csv('dn_ds_by_site_bovine_ha.tsv', sep='\t')

p = ggplot(mut_counts, aes(x=pos_aa, y=ratio)) +
  geom_point() +
  geom_text_repel(aes(label=ifelse(ratio>2.5,pos_aa,''))) +
  ggtitle('Selection in HA, bovine B3.13') +
  xlab('amino acid position') +
  ylab('dn/(ds+1)')

print(p)