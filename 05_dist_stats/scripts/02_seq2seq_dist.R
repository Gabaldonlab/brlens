library(e1071)
library(Hmisc)
library(corrplot)

library(ggplot2)
library(ggpubr)
theme_set(theme_bw())

stats <- function(x) {
  y <- c(summary(x, na.rm = TRUE), 'skew' = skewness(x, na.rm = TRUE),
         'kurt' = kurtosis(x, na.rm = TRUE), 'sum' = sum(x, na.rm = TRUE))
  return(y)
}

phy_files <- list.files(path = '../data/', pattern = '0005', full.names = TRUE)

remove(dat)
for (file in phy_files) {
  if (!exists('dat')) {
    dat <- read.csv(file)
  } else {
    dat <- rbind(dat, read.csv(file))
  }
}

write.csv(dat, '../outputs/0005_dists.csv')

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

write.csv(stats_df, file = '../outputs/0005_stats.csv')

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
  ylab('Duplications / speciations')

plot(dat[sample(which(dat$dist > 0 & dat$dist_norm < 20), 4000), 8:11],
     pch = 20)

dat_cor <- rcorr(as.matrix(dat[which(dat$dist > 0 & dat$dist_norm < 20), 8:11]))
corrplot(dat_cor$r, method = 'ellipse', type = 'lower', p.mat = dat_cor$P,
         sig.level = 0.05)
corrplot(add = TRUE, dat_cor$r, type = 'upper', method = 'number', tl.pos = 't')

stats_df <- data.frame(a[, , 1:4])
dist_df <- data.frame(sps, stats_df$Mean.dist)
spset <- dat$from_sp[!duplicated(dat$from_sp)]

b <- matrix(NA, nrow = length(spset), ncol = length(spset), dimnames = list(spset, spset))
for (i in 1:dim(dist_df)[1]) {
  b[dist_df[i, 1], dist_df[i, 2]] <- dist_df[i, 3]
  b[dist_df[i, 2], dist_df[i, 1]] <- dist_df[i, 3]
}

heatmap(b[-39, -39], na.rm = TRUE)
