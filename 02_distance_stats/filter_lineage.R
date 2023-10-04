library(tidyr)
library(stringr)

rdat <- read.csv('../../01_distances/outputs_lineage/mammal_dist_lineage.csv',
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

tdat <- read.csv('../../01_distances/sp_tree_lineage.csv',
                 check.names = FALSE)
ntdat <- tdat[, c(1, grep('_ndist', names(tdat)))]
ntdat <- gather(ntdat, 'node', 'ndist', -seed)
ntdat$node <- str_split(ntdat$node, '_', simplify = TRUE)[, 1]

rtdat <- tdat[, c(1, grep('_dist', names(tdat)))]
rtdat <- gather(rtdat, 'node', 'dist', -seed)
rtdat$node <- str_split(rtdat$node, '_', simplify = TRUE)[, 1]

rtdat <- merge(ntdat, rtdat, by = c('seed', 'node'))

mode <- function(x) {
  y <- density(x)
  y[['x']][which.max(y[['y']])]
}

y <- c(rtdat$dist)
x <- c(by(mgdat$ndist, mgdat$node, mode))

plot(x, y)

library(ggplot2)
theme_set(theme_bw())

ggplot() +
  geom_point(aes(x, y)) +
  geom_smooth(method = 'lm', data = data.frame(x, y), aes(x, y)) +
  xlab('Mode of the normalised distance') +
  ylab('Dated tree distance (Mya)')
