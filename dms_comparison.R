library(dplyr)
library(ggplot2)
library(ggrepel)

theme_set(theme_bw())
setwd('/Users/james/Documents/convergent/out')

mut_counts = read.csv('dn_ds_by_site_bovine_ha.tsv', sep='\t')

codons = read.csv('codon_mutation_counts_bovine_ha.tsv', sep='\t')

aminos = codons %>% 
  group_by(pos_aa, ref_aa, alt_aa) %>% 
  summarize()
  # mutate(impact = case_when(
  #   ref_aa == alt_aa ~ 'synonymous',
  #   .default = 'non-synonymous'
  # ))

dms = read.csv('../data/dms_data_ha.tsv', sep='\t')
dms = dms %>% 
  rename(
    pos_aa = position_aa,
    dms_name = name,
    dms_value = value,
  ) %>% 
  select(!c('ref_codon', 'alt_codon', 'gff_feature')) %>%
  group_by(ref_aa, alt_aa, pos_aa, dms_name, dms_value)

aminos_dms = left_join(aminos, dms, by=c('pos_aa', 'ref_aa', 'alt_aa'))

# what frac of amino mutations have dms values?
count_uq_pos = aminos_dms$pos_aa %>% unique() %>% length()
aminos_dms %>%
  ungroup() %>%
  select(c('dms_name', 'pos_aa')) %>%
  group_by(dms_name, pos_aa) %>% 
  group_by(dms_name) %>%
  summarize(n = n()) %>% 
  mutate(frac_positions = n / count_uq_pos) %>% 
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
  geom_point() +
  facet_wrap(~dms_name, scales='free', ncol=4) +
  geom_text_repel(aes(label=pos_aa), max.overlaps=5, size=4)
print(p)


