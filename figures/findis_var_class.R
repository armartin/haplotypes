## make figure 5A, figure S8, and figure S9
## many diagnostic plots towards end of script

library(dplyr)
library(cowplot)
library(RColorBrewer)
setwd('/Users/alicia/daly_lab/finrisk/findis')

pli <- read.table('fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt', header=T)

read_filter_data <- function(filename) {
  hap_exome <- read.table(gzfile(filename), header=T)
  hap_exome <- subset(hap_exome, pass=='true'&chrom %in% c(1:22))
  ## pLI > 0.95 for LoF, mis_z > 3.09 for missense
  hap_exome$freq <- (hap_exome$n_het  + 2* hap_exome$n_hom_alt) / (2*(hap_exome$n_hom_ref + hap_exome$n_het + hap_exome$n_hom_alt))
  hap_exome$ac <- hap_exome$n_het  + 2* hap_exome$n_hom_alt
  hap_exome$hh_ratio <- hap_exome$shared_het_1 / (hap_exome$shared_het_1 + hap_exome$non_shared_het)
  hap_exome$rr_ratio <- hap_exome$shared_ref_0 / (hap_exome$shared_ref_0 + hap_exome$non_shared_ref)
  hap_exome$hr_ratio <- hap_exome$hh_ratio / hap_exome$rr_ratio
  hap_exome$het_ref_enrich <- hap_exome$hh_ratio / hap_exome$rr_ratio
  hap_exome$var_class <- case_when(
    hap_exome$annotation %in% c('stop_gained', 'splice_acceptor_variant', 'splice_donor_variant') ~ 'LoF',
    hap_exome$annotation == 'synonymous_variant' ~ 'synonymous_variant',
    hap_exome$annotation == 'missense_variant' ~ 'missense_variant'
  )
  
  hap_exome$var_class_suffix <- case_when(
    hap_exome$annotation_suffix %in% c('stop_gained-HC', 'splice_acceptor_variant-HC', 'splice_donor_variant-HC') ~ 'LoF (high confidence)',
    hap_exome$annotation_suffix %in% c('stop_gained-LC', 'splice_acceptor_variant-LC', 'splice_donor_variant-LC') ~ 'LoF (low confidence)',
    hap_exome$annotation_suffix == 'synonymous_variant' ~ 'Synonymous',
    hap_exome$annotation_suffix == 'missense_variant-probably_damaging' ~ 'Missense (probably damaging)',
    hap_exome$annotation_suffix == 'missense_variant-possibly_damaging' ~ 'Missense (possibly damaging)',
    hap_exome$annotation_suffix == 'missense_variant-benign' ~ 'Missense (benign)',
    hap_exome$annotation_suffix == 'missense_variant' ~ 'Missense (unannotated)'
  )
  
  pli_hap <- merge(pli, hap_exome, by.x='gene', by.y='genes') %>%
    mutate(lof_constrained=pLI>0.95, missense_constrained=mis_z>3.09) %>%
    subset(!is.na(var_class_suffix)&var_class_suffix!='LoF (low confidence)')
  
  pli_hap$var_class_suffix <- factor(pli_hap$var_class_suffix, levels=c('LoF (high confidence)', 'Missense (probably damaging)', 'Missense (possibly damaging)',
                                                                        'Missense (benign)', 'Missense (unannotated)', 'Synonymous'))
  return(pli_hap)
}

pli_hap <- read_filter_data('daly_atgu_finnish_swedish_exomes.vep.haps.tsv.bgz')
hap_exome2 <- subset(pli_hap, !is.na(var_class_suffix))

brewer_vec <- c(brewer.pal(3, 'Reds')[3], brewer.pal(5, 'Blues')[5:2], brewer.pal(3, 'Set1')[3])
names(brewer_vec) <- sort(unique(hap_exome2$var_class_suffix))
  
p1 <- ggplot(hap_exome2, aes(x=freq, y=hh_ratio/rr_ratio, color=var_class_suffix, fill=var_class_suffix)) +
  geom_smooth(formula=y~log(x)) +
  scale_x_log10() +
  scale_color_manual(values=brewer_vec) +
  scale_fill_manual(values=brewer_vec) +
  labs(x='Allele frequency', y='Het:Ref pair sharing ratio') +
  theme_bw()# +
  
plot_overall <- function(dataset, title) {
  my_plot <- ggplot(dataset, aes(x=ac, y=hr_ratio, color=var_class_suffix, fill=var_class_suffix)) +
    geom_smooth(formula=y~log(x)) +
    coord_cartesian(ylim = c(-5, 23)) +
    scale_x_log10() +
    scale_x_log10(sec.axis = sec_axis(~./(9366*2), name = "Allele frequency", breaks=c(0.001, 0.01, 0.1, 1))) +
    scale_color_manual(values=brewer_vec, name='Variant class') +
    scale_fill_manual(values=brewer_vec, name='Variant class') +
    labs(x='Allele count', y='Het:Ref pair sharing ratio', title=title) +
    theme_bw() +
    theme(text = element_text(size=14),
          legend.position=c(0.4,0.3),
          legend.text=element_text(size=10),
          legend.background = element_rect(fill='transparent'))
  return(my_plot)
}
p2 <- plot_overall(subset(pli_hap, var_class_suffix != 'Missense (unannotated)'&ac<0.25*(9366*2)), 'All variants')
p2_nocpg <- plot_overall(subset(pli_hap, var_class_suffix != 'Missense (unannotated)'&cpg=='false'&ac<0.25*(9366*2)), 'No CpGs')
p2_nocpg2 <- p2 + labs(title=NULL) +
  theme(legend.position=c(0.3,0.2),
        text = element_text(size=18),
        legend.text = element_text(size=16),
        axis.text=element_text(colour = 'black'))
save(p2_nocpg2, file='het_ref_freq_nocpg.RData')
ggsave('het_ref_freq.pdf', p2 + labs(title=NULL), width=7, height=7)

plot_constraint <- function(dataset, title, constraint, name_constrain) {
  my_plot <- ggplot(dataset, aes_string(x='ac', y='hr_ratio', color='var_class_suffix', fill='var_class_suffix', linetype=constraint)) +
    geom_smooth(formula=y~log(x)) +
    labs(title=title) +
    coord_cartesian(ylim = c(-11, 25)) +
    #scale_x_log10(limits=c(1e-4,1e-1)) +
    scale_x_log10() +
    scale_x_log10(sec.axis = sec_axis(~./(9366*2), name = "Allele frequency", breaks=c(0.001, 0.01, 0.1, 1))) +
    scale_color_manual(values=brewer_vec, name='Variant class') +
    scale_fill_manual(values=brewer_vec, name='Variant class') +
    scale_linetype_manual(values=c('solid', 'dotted'), name=name_constrain) +
    guides(fill=F, color=F) +
    labs(x='Allele count', y='Het:Ref pair sharing ratio') +
    theme_bw() +
    theme(text = element_text(size=12),
          legend.position=c(0.3,0.15),
          legend.text=element_text(size=10),
          legend.background = element_rect(fill='transparent'))
}

p3_legend <- ggplot(subset(pli_hap, var_class_suffix!='Missense (unannotated)'), aes(x=ac, y=hh_ratio/rr_ratio, color=var_class_suffix, fill=var_class_suffix)) +
  geom_smooth(formula=y~log(x)) +
  scale_color_manual(values=brewer_vec, name='Variant class') +
  scale_fill_manual(values=brewer_vec, name='Variant class')

legend <- get_legend(p3_legend)# + theme(legend.position="bottom"))
pli_hap$hr_ratio <- pli_hap$hh_ratio / pli_hap$rr_ratio
p3 <- plot_constraint(subset(pli_hap, var_class_suffix!='Missense (unannotated)'), 'All variants, LoF constraint',
                      'lof_constrained', 'LoF constrained')
p3_nocpg <- plot_constraint(subset(pli_hap, var_class_suffix!='Missense (unannotated)'&cpg=='false'), 
                            'No CpGs, LoF constraint', 'lof_constrained', 'LoF constrained')
p4 <- plot_constraint(subset(pli_hap, !var_class_suffix %in% c('LoF (high confidence)', 'Missense (unannotated)')), 
                      'All variants, Missense constraint', 'missense_constrained', 'Missense constrained')
p4_nocpg <- plot_constraint(subset(pli_hap, !var_class_suffix %in% c('LoF (high confidence)', 'Missense (unannotated)')&cpg=='false'),
                            'No CpGs, Missense constraint', 'missense_constrained', 'Missense constrained')

p3_p4_plot <- plot_grid(p3, p3_nocpg, NULL, p4, p4_nocpg, NULL, ncol=3, labels=c('A', 'B', '', 'C', 'D', ''), 
                        rel_widths=c(1,1,0.75))
p3_p4_plot <- p3_p4_plot + draw_grob(legend, x=2.1/2.75, y=0.25, width=0.2, height=0.5)
save_plot('het_ref_freq_cpg_constraint_comparison.pdf', p3_p4_plot, base_height=8, base_width=10)


# Plot by haplotype length ------------------------------------------------
setwd('/Users/alicia/daly_lab/finrisk/haplotypes/eagle/')

pli_hap1 <- read_filter_data('daly_atgu_finnish_swedish_exomes.vep.haps.1.0-1.5cM.tsv.bgz')
pli_hap2 <- read_filter_data('daly_atgu_finnish_swedish_exomes.vep.haps.1.5-2.0cM.tsv.bgz')
pli_hap3 <- read_filter_data('daly_atgu_finnish_swedish_exomes.vep.haps.2.0-2.5cM.tsv.bgz')
pli_hap4 <- read_filter_data('daly_atgu_finnish_swedish_exomes.vep.haps.2.5-3.0cM.tsv.bgz')
pli_hap5 <- read_filter_data('daly_atgu_finnish_swedish_exomes.vep.haps.3.0-3.5cM.tsv.bgz')
pli_hap6 <- read_filter_data('daly_atgu_finnish_swedish_exomes.vep.haps.3.5-4.0cM.tsv.bgz')
pli_hap7 <- read_filter_data('daly_atgu_finnish_swedish_exomes.vep.haps.4.0-4.5cM.tsv.bgz')
pli_hap8 <- read_filter_data('daly_atgu_finnish_swedish_exomes.vep.haps.4.5-5.0cM.tsv.bgz')

pli_hap1$interval <- '1.0-1.5 cM'
pli_hap2$interval <- '1.5-2.0 cM'
pli_hap3$interval <- '2.0-2.5 cM'
pli_hap4$interval <- '2.5-3.0 cM'
pli_hap5$interval <- '3.0-3.5 cM'
pli_hap6$interval <- '3.5-4.0 cM'
pli_hap7$interval <- '4.0-4.5 cM'
pli_hap8$interval <- '4.5-5.0 cM'

pli_hap_intervals <- bind_rows(pli_hap1, pli_hap2, pli_hap3, pli_hap4, pli_hap5, pli_hap6, pli_hap7, pli_hap8)

my_plot <- ggplot(subset(pli_hap_intervals, var_class_suffix != 'Missense (unannotated)'&ac<0.25*(9366*2)), aes(x=ac, y=hr_ratio, color=var_class_suffix, fill=var_class_suffix)) +
  facet_wrap(~interval, ncol=4) +
  geom_smooth(formula=y~log(x)) +
  scale_x_log10() +
  scale_color_manual(values=brewer_vec, name='Variant class') +
  scale_fill_manual(values=brewer_vec, name='Variant class') +
  labs(x='Allele count', y='Het:Ref pair sharing ratio') +
  theme_bw() +
  theme(text = element_text(size=14),
        legend.text=element_text(size=10),
        legend.background = element_rect(fill='transparent'))

ggsave('het_ref_hap_intervals.pdf', my_plot, width=10, height=5)

pli_hap2 <- pli_hap %>%
  subset(freq <=0.001) %>%
  mutate(ac=n_het+2*n_hom_alt) %>%
  group_by(ac, var_class_suffix) %>%
  summarize(mean=mean(hh_ratio/rr_ratio, na.rm=T), sd=sd(hh_ratio/rr_ratio, na.rm=T), freq=mean(freq))

p3_test <- ggplot(subset(pli_hap, freq>=0.001), aes(x=freq, y=hh_ratio/rr_ratio, color=var_class_suffix, fill=var_class_suffix)) +#, linetype=lof_constrained)) +
  geom_smooth(formula=y~log(x)) +
  geom_point(data=pli_hap2, aes(x=freq, y=mean)) +
  labs(title='CpGs included, LoF constraint') +
  scale_x_log10() +
  scale_color_manual(values=brewer_vec) +
  scale_fill_manual(values=brewer_vec) +
  labs(x='Allele frequency', y='Het:Ref pair sharing ratio') +
  theme_bw()# +

pli_hap %$% table(cpg, var_class_suffix, ac==2)
p3_test2 <- ggplot(subset(pli_hap, cpg=='false'), aes(x=ac, y=hh_ratio/rr_ratio, color=var_class_suffix, fill=var_class_suffix)) +#, linetype=lof_constrained)) +
  geom_smooth(formula=y~log(x)) +
  scale_x_log10(sec.axis = sec_axis(~./(9366*2), name = "Allele frequency", breaks=c(0.001, 0.01, 0.1, 1))) +
  labs(title='CpGs excluded, LoF constraint') +
  scale_x_log10() +
  scale_color_manual(values=brewer_vec) +
  scale_fill_manual(values=brewer_vec) +
  labs(x='Allele count', y='Het:Ref pair sharing ratio') +
  theme_bw()# +
ggsave('p3_2.pdf', p3_test2)

my_plot <- plot_grid(p1, p2, p3, p4, labels=c('A', 'B', 'C', 'D'))
save_plot('het_ref_freq_nocpg.pdf', my_plot, nrow=2, base_height=6, base_width=14)

p2_plot <- plot_grid(p2, p2_nocpg, labels=c('A', 'B'), align='h')
save_plot('het_ref_freq_cpg_comparison.pdf', p2_plot, nrow=2, base_height=2, base_width=7)

ggsave('carrier_ref_freq.pdf', p1, height=5, width=7)
ggsave('carrier_ref_freq_class.pdf', p2, height=5, width=7)
findis <- read.table(gzfile('findis_merge2.summary.txt.gz'), header=T)
findis$freq <- (findis$a_tot*2 + findis$h_tot) / (2*(findis$r_tot + findis$h_tot + findis$a_tot))
findis$rr_ratio <- findis$rr_haps / findis$rr_tot
findis$hh_ratio <- findis$hh_haps / findis$hh_tot
findis$h_r_ratio <- findis$hh_ratio / findis$rr_ratio
findis2 <- subset(findis, freq < 0.05 & freq > 0)

subset(findis2, pos %in% c(178359918, 91449319, 107427289, 15538697))
ggplot(findis2, aes(x=freq, y=hh_ratio)) +
  geom_point() +
  geom_smooth()
  
ggplot(findis2, aes(x=freq, y=h_r_ratio)) +
  geom_point() +
  geom_smooth(method='lm')

ggplot(hap_exome_compare, aes(x=freq, y=het_ref_enrich, color=var_class)) +
  geom_smooth(formula=y~log(x))

ggplot(hap_exome_compare, aes(x=freq, y=hh_ratio, color=var_class)) +
  geom_smooth(formula=y~log(x)) +
  labs(x='Allele frequency', y='Carrier sharing ratio')

ggplot(hap_exome_compare, aes(x=freq, y=rr_ratio, color=var_class)) +
  geom_smooth(formula=y~log(x)) +
  labs(x='Allele frequency', y='Ref sharing ratio')

p1 <- ggplot(hap_exome_compare, aes(x=freq, y=hh_ratio/rr_ratio)) +
  geom_smooth(formula=y~log(x)) +
  labs(x='Allele frequency', y='Carrier:Ref sharing ratio')
