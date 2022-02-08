library(e1071)
library(ggplot2)
library(ggpubr)
theme_set(theme_bw())

stats <- function(x) {
  y <- c(summary(x, na.rm = TRUE), 'skew' = skewness(x, na.rm = TRUE),
         'kurt' = kurtosis(x, na.rm = TRUE), 'sum' = sum(x, na.rm = TRUE))
  return(y)
}
dat <- read.csv('../data/0435_dist.csv')


sps <- t(combn(dat$from_sp[!duplicated(dat$from_sp)], 2))
combs <- apply(sps[, 1:2], 1, paste, collapse = ' - ')

choose(length(dat$from_sp[!duplicated(dat$from_sp)]), 2)

remove(a)
# pdf('~/Documents/pairs_descr.pdf', width = 10, height = 5)
for (i in 1:dim(sps)[1]) {
  print(i)
  if (sum(dat$from_sp == sps[i, 1] & dat$to_sp == sps[i, 2]) > 0) {
    sdat <- dat[which(dat$from_sp == sps[i, 1] & dat$to_sp == sps[i, 2]), ]
    sdat_stats <- apply(sdat[, 8:11], 2, stats)
    if (!exists('a')) {
      a <- array(dim = c(dim(sps)[1], dim(sdat_stats)[1], 4),
                 dimnames = list(combs, names(stats(rnorm(100))),
                                 c('dist', 'dist_norm', 'sp', 'dupl')))
    }
    a[i, , 'dist'] <- sdat_stats[, 1]
    a[i, , 'dist_norm'] <- sdat_stats[, 2]
    a[i, , 'sp'] <- sdat_stats[, 3]
    a[i, , 'dupl'] <- sdat_stats[, 4]
  }
  
  # pa <- ggplot(sdat, aes(dist)) +
  #   geom_density() +
  #   xlab('Distance') +
  #   labs(title = combs[i])
  # 
  # pb <- ggplot(sdat, aes(sp)) +
  #   geom_histogram() +
  #   xlab('Speciation events') +
  #   labs(title = combs[i])
  # 
  # pc <- ggplot(sdat, aes(dupl)) +
  #   geom_histogram() +
  #   xlab('Duplication events') +
  #   labs(title = combs[i])
  # 
  # pd <- ggplot(sdat, aes(dist, factor(dupl))) +
  #   geom_violin(draw_quantiles = c(0.5)) +
  #   xlab('Distance') +
  #   ylab('Duplications') +
  #   labs(title = combs[i])
  # 
  # pe <- ggplot(sdat, aes(dist, factor(sp))) +
  #   geom_violin(stat = 'ydensity', draw_quantiles = c(0.5)) +
  #   xlab('Distance') +
  #   ylab('Speciations') +
  #   labs(title = combs[i])
  # 
  # pf <- ggplot(sdat, aes(dist, dupl / sp)) +
  #   geom_point() +
  #   xlab('Distance') +
  #   ylab('Duplications / speciations') +
  #   labs(title = combs[i])
  # 
  # print(ggarrange(pa, pb, pc, pd, pe, pf, align = 'hv'))
}
# dev.off()

stats_df <- data.frame(a[, , 1:4])
stats_df <- na.omit(stats_df)

str(stats_df)

ggplot(stats_df, aes(Median.dist)) +
  geom_density()

ggplot(stats_df, aes(Median.dist_norm)) +
  geom_density()

ggplot(stats_df, aes(Mean.dist)) +
  geom_density()

ggplot(stats_df, aes(Mean.dist_norm)) +
  geom_density()

ggplot(stats_df, aes(Median.dupl)) +
  geom_density()

ggplot(stats_df, aes(Mean.dupl)) +
  geom_density()

ggplot(stats_df, aes(Median.sp)) +
  geom_density()

ggplot(stats_df, aes(Mean.sp)) +
  geom_density()

ggplot(stats_df, aes(skew.dist)) +
  geom_density()

ggplot(stats_df, aes(Mean.dist_norm, Mean.sp)) +
  geom_point()

ggplot(stats_df, aes(Median.dist_norm, Median.sp)) +
  geom_point()

ggplot(stats_df, aes(Mean.dist, Mean.dupl)) +
  geom_point() +
  xlab('Mean distance') +
  ylab('Mean number duplications')

ggplot(stats_df, aes(Mean.dist_norm, Mean.dupl)) +
  geom_point() +
  xlab('Mean nrom distance') +
  ylab('Mean number duplications')

ggplot(stats_df, aes(Median.dist, Median.dupl)) +
  geom_point()

ggplot(stats_df, aes(Median.dist_norm, sum.dupl / sum.sp)) +
  geom_point() +
  xlab('Distance') +
  ylab('Duplications / speciations') +
  labs(title = combs[i])
