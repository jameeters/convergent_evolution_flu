library(dplyr)
library(ggplot2)
library(ggrepel)

theme_set(theme_bw())
setwd('/Users/james/Documents/convergent')

mut_counts = read.csv('out/dn_ds_by_site_bovine_ha.tsv', sep='\t') %>% 
  mutate(host='bovine')

mut_counts = full_join(
  mut_counts,
  read.csv('out_avian_b3.13/dn_ds_by_site_ha.tsv', sep='\t') %>% 
    mutate(host='avian')
)

mut_counts = mut_counts %>% 
  mutate(
    ratio = case_when(
      count_s == 0 ~ count_n,
      .default = count_n / count_s
    ),
    zero_synonymous = (count_s == 0)
  )

codons = read.csv('out/codon_mutation_counts_bovine_ha.tsv', sep='\t') %>% 
  mutate(host='bovine')
codons = full_join(
  codons,
  read.csv('out_avian_b3.13/codon_mutation_counts_ha.tsv', sep='\t') %>% 
    mutate(host='avian')
)

aminos = codons %>% 
  group_by(pos_aa, ref_aa, alt_aa, host) %>% 
  summarize()

dms = read.csv('data/dms_data_ha.tsv', sep='\t')
dms = dms %>% 
  rename(
    pos_aa = position_aa,
    dms_name = name,
    dms_value = value,
  ) %>% 
  select(!c('ref_codon', 'alt_codon', 'gff_feature')) %>%
  group_by(ref_aa, alt_aa, pos_aa, dms_name, dms_value)

aminos_dms = left_join(aminos, dms, by=c('pos_aa', 'ref_aa', 'alt_aa'))

# *********************************************
# what frac of amino mutations have dms values?
# *********************************************
count_uq_pos = aminos_dms %>% 
  distinct(pos_aa, host) %>% 
  group_by(host) %>% 
  summarize(n_uq_pos=n())
  
aminos_dms %>%
  ungroup() %>%
  select(c('dms_name', 'pos_aa', 'host')) %>%
  group_by(dms_name, pos_aa, host) %>%
  group_by(dms_name, host) %>%
  summarize(n = n()) %>%
  left_join(count_uq_pos, by='host') %>% 
  mutate(frac_positions = n / n_uq_pos ) %>%
  print()

aminos_dms = aminos_dms %>% 
  group_by(pos_aa, dms_name) %>% 
  summarize(
    mean_dms_value = mean(dms_value),
    med_dms_value = median(dms_value),
    min_dms_value = min(dms_value),
    max_dms_value = max(dms_value),
  )

dn_ds_dms = mut_counts %>%
  left_join(aminos_dms, by=c('pos_aa')) %>% 
  filter(!is.na(dms_name)) %>% 
  filter(!dms_name %in% c('evescape_sigmoid', 'evescape'))

p = ggplot(dn_ds_dms, aes(x=mean_dms_value, y=ratio)) +
  geom_point(aes(shape=zero_synonymous)) +
  ggtitle('HA, B3.13') +
  ylab('dn/ds') +
  facet_wrap(vars(dms_name, host), scales='free_x', ncol=4) +
  geom_text_repel(aes(label=pos_aa), max.overlaps=5, size=4)
print(p)


