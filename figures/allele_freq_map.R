#!/usr/bin/env Rscript

library(geosphere)
library(ggmap)
library(ggplot2)
library(dplyr)
library(rgdal)
library(RColorBrewer)
library(mapdata)
library(fields)
library(data.table)
library(argparse)

## get arguments
parser <- ArgumentParser()
parser$add_argument('--disease', type='character')
parser$add_argument('--variant', type='character')
args <- parser$parse_args()

print(paste(args$disease, args$variant, sep=', '))
## load exomes and birth record data
setwd('/Users/alicia/daly_lab/finrisk/findis')
exomes <- read.table(gzfile('findis_genotypes.tsv.gz'), header=T, sep='\t')
#diar1 <- exomes[grepl('7:107427289', exomes$VARIANT),]
disease <- exomes[grepl(args$variant, exomes$VARIANT),]
freq_table <- table(disease$GENO)
table_names <- as.numeric(names(table(disease$GENO)))
freq <- weighted.mean(freq_table, table_names) / (2 * sum(freq_table))
print(paste(args$disease, args$variant, freq, sep=', '))
birth <- read.csv('../data/all_birth_wgs_exome_combined.csv', header=T)
exome_disease <- merge(birth, disease, by.x='exome_id', by.y='IID')
exome_disease <- subset(exome_disease, Birth.records.avail==1)

## load map and municipality info
load('../data/maps/suomi.df.RData')

## find scenarios where parents born within 50 miles of each other, store mean parent locations
mun_data <- exome_disease %>%
  subset(!is.na(KUNT_K14)) %>%
  left_join(mun_locs, by=c('KUNT_K14' = 'CODE')) %>%
  select(-c(contains('NAME'), contains('DETAILS'), ends_with('LAAN'))) %>%
  rename(LAT_mom=LAT, LON_mom=LON) %>%
  left_join(mun_locs, by=c('KUNT_K15' = 'CODE')) %>%
  select(-c(contains('NAME'), contains('DETAILS'), ends_with('LAAN'))) %>%
  rename(LAT_dad=LAT, LON_dad=LON) %>%
  subset(!is.na(LAT_dad) & !is.na(LAT_mom)) %>%
  mutate(parent_dist = distCosine(cbind(LON_mom, LAT_mom), cbind(LON_dad, LAT_dad))) %>%
  subset(parent_dist < 80000) %>%
  mutate(meanp_lat = (LAT_mom + LAT_dad)/2, meanp_lon = (LON_mom + LON_dad)/2)

pair_dist <- rdist.earth(select(mun_data, meanp_lon, meanp_lat), select(mun_locs, LON, LAT))
mun_data$mun <- mun_locs$NAME[(apply(pair_dist, 1, which.min))]
mun_data$meanp_code <- mun_locs$CODE[(apply(pair_dist, 1, which.min))]

## make regional data
reg_data <- subset(exome_disease, FR_cohort == 'FR07')
reg_data <- reg_data %>%
  full_join(mun_locs, by=c('SKUNTA_LAAN'='LAAN')) %>%
  mutate(mun=NAME)
  
reg_data <- reg_data %>%
  group_by(ID2) %>% 
  mutate(weight=1/n()) %>%
  ungroup()

mun_data$weight <- 1

freq_summary <- bind_rows('mun'=mun_data, 'reg'=reg_data, .id='group') %>%
  #group_by(max_lat, max_lon, min_lat, min_lon) %>%
  group_by(mun) %>%
  dplyr::summarize(freq = weighted.mean(GENO, weight), n=n()) 

freq_summary <- exome_disease %>%
  group_by(SKUNTA_LAAN) %>%
  dplyr::summarize(freq = mean(GENO/2))

#suomi.df$freq <- freq_summary$freq[match(as.character(suomi.df$id), as.character(freq_summary$mun))]
suomi.df$freq <- freq_summary$freq[match(as.character(suomi.df$mun), as.character(freq_summary$SKUNTA_LAAN))]
#suomi.df$n <- freq_summary$n[match(as.character(suomi.df$id), as.character(freq_summary$n))]

suomi.df$disease <- args$disease

save(suomi.df, file=paste0(args$disease, '_', args$variant, '.RData'))