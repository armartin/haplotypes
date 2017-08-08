#!/usr/bin/env Rscript

library(argparse)
library(RColorBrewer)
devtools::install_github('konradjk/MapGAM')
library(MapGAM)
library(mapplots)
library(geosphere)
library(plyr)
library(dplyr)

parser <- ArgumentParser()
parser$add_argument('--point_num', type='integer')
args <- parser$parse_args()

setwd('/home/unix/armartin/finrisk/haplotypes/maps')
grid_points <- read.table('fin_grid.txt', header=T)
my_point = grid_points[args$point_num,] #c(28.975840, 65.924915) #random place
## load ibd data
cum_ibd <- read.table(gzfile('/home/unix/armartin/finrisk/haplotypes/ibd/Engagex_Egfext_Migraine_FinnishCardio_FTC.cum_pairs.ibd.gz'), header=T)
cum_ibd <- subset(cum_ibd, !is.na(lat1) & !is.na(lat2))
cum_ibd$dist <- distCosine(cbind(cum_ibd$lon1, cum_ibd$lat1), cbind(cum_ibd$lon2, cum_ibd$lat2)) #meters
# cum_ibd_close <- subset(cum_ibd, dist < 80e3)
# cum_ibd_close$mean_lat <- apply(cbind(cum_ibd_close$lat1, cum_ibd_close$lat2), 1, mean)
# cum_ibd_close$mean_lon <- apply(cbind(cum_ibd_close$lon1, cum_ibd_close$lon2), 1, mean)
# cum_ibd_close$cum_ibd_mb <- cum_ibd_close$cum_ibd / 1e6
# cum_ibd_close <- subset(cum_ibd_close, cum_ibd_mb < 350)
set.seed(42)
#cum_ibd_thin <- sample_n(cum_ibd_close, 10000, replace=F)

## make a grid of finland
# xmin=min(grid_points$max_lon)
# xmax=max(grid_points$max_lon)
# ymin=min(grid_points$max_lat)
# ymax=max(grid_points$max_lat)
# grd <- data.frame(make.grid(cum_ibd_thin$lon1, cum_ibd_thin$lat1, cum_ibd_thin$cum_ibd_mb, byx = 0.2, byy = 0.1, 
#                             xlim = c(xmin-0.7, xmax+0.7), ylim = c(ymin, ymax+1.5)))
# draw.grid(grd)
# grd$val <- rownames(grd)
# long_grid <- melt(grd, id="val")
# colnames(long_grid) <- c('max_lon', 'max_lat', 'useless')
# long_grid$max_lat <- as.numeric(gsub('X', '', long_grid$max_lat))
# long_grid$max_lon <- as.numeric(long_grid$max_lon)

## get distance between each individual and grid point.
#most_ibd <- map_fit$grid[which(map_fit$fit==max(map_fit$fit)),]
cum_ibd$dist_highest1 <- distCosine(cbind(cum_ibd$lon1, cum_ibd$lat1), my_point)
cum_ibd$dist_highest2 <- distCosine(cbind(cum_ibd$lon2, cum_ibd$lat2), my_point)
cum_ibd_highest <- subset(cum_ibd, dist_highest1 < 50000 | dist_highest2 < 50000) #get only pairs where at least 1 individual is within 80 km
#helsinki: 1,348,177
#random: 158,617
#igloo: 3985
cum_ibd_highest$which_max <- as.numeric(apply(as.matrix(cum_ibd_highest[,c("dist_highest1","dist_highest2")]), 1, which.max))
cum_ibd_highest <- adply(cum_ibd_highest, 1, function(x)  #calculates max position between pairs. naming is mean (inconsistent) to be consistent with my_grid
{data.frame(max_lat=as.numeric((x[2*as.numeric(x['which_max'])+2] )),
            max_lon=as.numeric((x[2*as.numeric(x['which_max'])+3] )),
            cum_ibd_mb=x['cum_ibd']/1e6)
})
#igloo: 3.669
#random: 154.392
cum_ibd_highest$cum_ibd_mb <- cum_ibd_highest$cum_ibd/1e6
cum_ibd_highest_thin <- sample_n(cum_ibd_highest, 10000, replace=F)
#cum_ibd_highest_thin <- cum_ibd_highest

# load map of Finland
load(url("http://biogeo.ucdavis.edu/data/gadm2/R/FIN_adm0.RData"))
suomi <- gadm
my_grid <- predgrid(grid_points, suomi)

# Predict map values in a grid and plot
date(); high_ibd_fit <- modgam(cum_ibd_highest[,c('cum_ibd_mb', 'max_lon', 'max_lat')], my_grid, m='crude', family='gaussian', sp=0.5); date()
#time to fit without thinning (random): ~15 minutes, angry laptop, R throwing errors.
#date(); map_fit <- modgam(cum_ibd_thin[,c(11,10,9)], my_grid, m='crude', family='gaussian', sep=0.2); date()

# generate fitted map of finland
png(paste0('point', args$point_num, '.png'), height=837, width=450, res=90, type='cairo')
colormap(high_ibd_fit, suomi, axes=T, col=colorRampPalette(rev(brewer.pal(9, 'YlOrBr'))), leglab='Cumulative IBD (Mb)', xlab='test', mar=c(8,5,5,3))#, mapmin=28.4, mapmax=118) #map needs to be produced from maps or maptools
points(my_point, pch=4, cex=3)
mtext(paste0("Point ", args$point_num), 3, line=2.5, font=2, cex=1.5)
dev.off()


# date(); png('point1.png', height=857, width=450, res=90)
# colormap(high_ibd_fit, suomi, axes=T, col=colorRampPalette(rev(brewer.pal(9, 'YlOrBr'))), leglab='Cumulative IBD (Mb)', xlab='test', mar=c(12,5,5,3)) #map needs to be produced from maps or maptools
# mtext("Close birth regions", 3, line=2.5, font=2, cex=1.5)
# dev.off(); date()
