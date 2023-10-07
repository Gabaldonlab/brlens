library(tidyr)
library(stringr)

get_spto <- function(x) {
  x[which(x != 'HUMAN')]
}

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

edat <- mgdat

save(edat, file = 'outputs/ev_dat.RData')

rdat <- read.csv('../01_distance_calculations/outputs/lng_dist.csv',
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

ldat <- mgdat

save(ldat, file = 'outputs/lng_dat.RData')

rdat <- read.csv('../01_distance_calculations/outputs/t2t_dist.csv')

datc <- rdat[which(rdat$MRCA_type == 'S' &
                     (rdat$from_seq == 'HUMAN' | rdat$to_seq == 'HUMAN') &
                     (rdat$from == rdat$tree | rdat$to == rdat$tree)), ]

datc[, 'sp_to'] <- apply(datc[, c(3, 5)], 1, get_spto)

tdat <- datc[which(datc$ndist <= quantile(datc$ndist, 0.9)), ]
tdat <- tdat[, c('tree', 'sp_to', 'dist', 'ndist')]
names(tdat) <- c('seed', 'node', 'dist', 'ndist')

save(tdat, file = 'outputs/t2t_dat.RData')

