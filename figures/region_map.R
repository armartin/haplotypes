## plot figure S1

library(geosphere)
library(ggmap)
library(ggplot2)
library(dplyr)
library(rgdal)
library(RColorBrewer)
library(mapdata)
library(fields)
library(data.table)

## load exomes and birth record data
setwd('/Users/alicia/daly_lab/finrisk/data/')

## load map and municipality info
load('maps/suomi.df.RData')

suomi.df.centroids <- data.frame(long=coordinates(suomi.spdf[,1]),
                                 lat=coordinates(suomi.spdf[,2]),
                                 MUN=suomi.spdf[,'MUN'])

centroids <- suomi.df.centroids %>%
  group_by(MUN) %>%
  mutate(mean_lat=mean(lat.2), mean_long=mean(long.1)) %>%
  ungroup() %>%
  select(MUN, mean_lat, mean_long) %>%
  #ungroup() %>%
  distinct() %>%
  arrange(MUN) %>%
  subset(MUN < 11)
# adjust centroid positions for 2 regions
add_centroids <- data.frame(MUN=c(11,12), mean_lat=c(64.75,67.4), mean_long=c(27,26.5046))
centroids <- bind_rows(centroids, add_centroids)

my_map <- ggplot(suomi.df, aes(x=long, y=lat, group=group)) +
  geom_polygon(aes(color=type, size=type, alpha=level), fill=brewer.pal(3, 'Set1')[2]) +
  scale_size_manual(values=size_array) +
  scale_color_manual(values=color_array) +
  xlab('Longitude') +
  ylab('Latitude') +
  geom_text(data=centroids, aes(x=mean_long, y=mean_lat, label=MUN, group=MUN), color='white') +
  labs(fill='Allele\nfrequency') +
  guides(alpha=F, color=F, size=F) +
  theme_bw() +
  theme(text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5))

date(); ggsave('mun_reg.png',
               my_map,
               width=2.25,
               height=3.5,
               units='in',
               dpi=300); date()

color <- c(rep('black', 5), rep('white', 6))
names(color) <- c(1,2,4:12)
my_map2 <- ggplot(suomi.df, aes(x=long, y=lat, group=group)) +
  geom_polygon(aes(color=type, size=type, alpha=level, fill=as.character(mun))) +
  scale_size_manual(values=size_array) +
  scale_color_manual(values=color_array) +
  scale_fill_manual(values=color_array) +
  xlab('Longitude') +
  ylab('Latitude') +
  geom_text(data=centroids, aes(x=mean_long, y=mean_lat, label=MUN, group=MUN), color=color) +
  labs(fill='Allele\nfrequency') +
  guides(alpha=F, color=F, size=F, fill=F) +
  theme_bw() +
  theme(text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5))

date(); ggsave('supp_mun_reg.png',
               my_map2,
               width=2.25,
               height=3.5,
               units='in',
               dpi=300); date()
