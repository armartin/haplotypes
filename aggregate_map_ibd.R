library(plyr)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
setwd('/Users/alicia/daly_lab/finrisk/data/maps')

cities <- data.frame(city=c('Helsinki', 'Turku', 'Tampere', 'Vaasa', 'Oulu', 'Rovaniemi', 'Kuopio', 'Ilomantsi', 'Kuusamo'),
                      region=c(rep('Major cities', 3), rep('Western cities', 3), rep('Eastern cities', 3)))

read_shapes <- function(city) {
  #print(city)
  suomi.df <- get(load(paste0('cumulative_ibd_', city, '.RData')))
  #suomi.df$mean_cum_ibd <- suomi.df$mean_cum_ibd[[1]] #this was messed up in generation, so fix
  #suomi.df$n_pairs <- suomi.df$n_pairs[[1]]
  #print(colnames(city_df))
  #print(tbl_df(city_df))
  return(suomi.df)
}

#test <- read_shapes('Tampere')
#test$mean_cum_ibd[[1]]

all.df <- ldply(cities$city, read_shapes)
all.df <- all.df %>%
  left_join(cities, by=c('loc_name'='city'))
all.df$loc_name <- factor(all.df$loc_name, levels=cities$city)

size_array <- c('region'=0.5, 'municipality'=0.1)
color_array1 <- colorRampPalette(brewer.pal(9, 'YlOrBr'))(12)
names(color_array1) <- 1:12
color_array <- c('region'='black', 'municipality'='grey20', 'suomi'='darkred', 'mun'='green',
                 color_array1)
my_points <- read.table('~/daly_lab/finrisk/data/finnish_cities.txt', header=T)
my_points <- my_points %>% left_join(cities, by=c('City'='city')) %>%
  mutate(loc_name=City)

plot_map <- function(which_region, xplot, yplot) {
  print(date())
  facet_maps <- ggplot(subset(all.df, region==which_region), aes(x=long, y=lat, group=group)) +
    geom_polygon(aes(color=type, size=type, fill=mean_cum_ibd, alpha=level)) +
    facet_grid(region ~ loc_name, scales='free') +
    scale_fill_gradientn(colours=rev(brewer.pal(11, 'RdYlBu'))) +#, limits=c(0,80)) + #set limits
    scale_size_manual(values=size_array) +
    scale_color_manual(values=color_array) +
    #ylab('Latitude') +
    #ggtitle(args$location) +
    labs(fill='Mean\npairwise\ncumulative\nIBD') +
    geom_point(data=subset(my_points, region==which_region), aes(x=Longitude, y=Latitude), group='2', shape='*', size=12, color=brewer.pal(3, 'Oranges')[3]) +
    guides(alpha=F, color=F, size=F) +
    theme_bw() +
    theme(text = element_text(size = 16),
          axis.text = element_text(size=10),
          plot.title = element_text(hjust = 0.5),
          legend.title = element_text(size=12),
          legend.text = element_text(size=10))
  if(xplot) {
    facet_maps = facet_maps +
      xlab('Longitude')
  } else {
    facet_maps = facet_maps + 
      xlab('') +
      theme(axis.text.x=element_blank())
  }
  if(yplot) {
    facet_maps = facet_maps + 
      ylab('Latitude')
  } else {
    facet_maps = facet_maps + 
      ylab('')
  }
  return(facet_maps) 
}

p1 <- plot_map('Major cities', FALSE, TRUE)
p2 <- plot_map('Western cities', FALSE, TRUE)
p3 <- plot_map('Eastern cities', TRUE, TRUE)

p_aggregate <- plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), ncol=1, align='v')
date(); ggsave(paste0('/Users/alicia/daly_lab/finrisk/data/maps/cumulative_ibd_facet_3cM_v3.png'),
       p_aggregate,
       width=6.5,
       height=9,
       units='in',
       dpi=300); date()
