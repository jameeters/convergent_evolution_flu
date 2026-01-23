library(dplyr)
library(ggplot2)
library(ggrepel)

theme_set(theme_bw())
setwd('/Users/james/Documents/convergent')

mut_counts_bovine = read.csv('out/dn_ds_by_site_bovine_ha.tsv', sep='\t')
mut_counts_avian = read.csv('out_avian_b3.13/dn_ds_by_site_ha.tsv', sep='\t')

mut_counts = full_join(
  mut_counts_bovine %>% mutate(host='bovine'),
  mut_counts_avian %>% mutate(host='avian')
)

mut_counts = mut_counts %>% 
  mutate(
    ratio = case_when(
      count_s == 0 ~ count_n,
      .default = count_n / count_s
    ),
    zero_synonymous = (count_s == 0)
  )

p = ggplot(mut_counts, aes(x=pos_aa, y=ratio)) +
  geom_point(aes(color=zero_synonymous)) +
  facet_wrap(~host, ncol=1) +
  geom_text_repel(aes(label=ifelse(ratio>2.5,pos_aa,''))) +
  ggtitle('Selection in HA, B3.13') +
  xlab('amino acid position') +
  ylab('dn/ds')
print(p)

