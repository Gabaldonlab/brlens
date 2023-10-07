library(tidyr)
library(stringr)

rdat <- read.csv('../01_distance_calculations/outputs/evt_dist.csv',
                 check.names = FALSE)

ncdat <- rdat[, c(1, grep('_ndist', names(rdat)))]
ngdat <- gather(ncdat, 'node', 'ndist', -seed)
ngdat$node <- str_split(ngdat$node, '_', simplify = TRUE)[, 1]

rcdat <- rdat[, c(1, grep('_dist', names(rdat)))]
rgdat <- gather(rcdat, 'node', 'dist', -seed)
rgdat$node <- str_split(rgdat$node, '_', simplify = TRUE)[, 1]

mgdat <- merge(rgdat, ngdat, by = c('seed', 'node'))

mgdat <- na.omit(mgdat)
mgdat <- mgdat[which(mgdat$ndist <= quantile(mgdat$ndist, 0.99)), ]

boxplot(mgdat$ndist ~ mgdat$node)
ldat <- mgdat

table(ldat$node)

save(ldat, file = '../data/lineage_dat.RData')
