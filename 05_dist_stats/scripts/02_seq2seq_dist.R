library(e1071)
library(ggplot2)
theme_set(theme_bw())

stats <- function(x) {
  y <- c(summary(x, na.rm = TRUE), 'skew' = skewness(x, na.rm = TRUE),
         'kurt' = kurtosis(x, na.rm = TRUE))
  return(y)
}

dat <- read.csv('../data/0003_dist.csv')


sps <- t(combn(dat$from_sp[!duplicated(dat$from_sp)], 2))
combs <- apply(sps[, 1:2], 1, paste, collapse = ' - ')

for (i in 1:dim(sps)[1]) {
  print(i)
  if (sum(dat$from_sp == sps[i, 1] & dat$to_sp == sps[i, 2]) > 0) {
    sdat <- dat[which(dat$from_sp == sps[i, 1] & dat$to_sp == sps[i, 2]), ]
    sdat_stats <- apply(sdat[, 8:10], 2, stats)
    if (!exists('a')) {
      a <- array(dim = c(dim(sps)[1], dim(sdat_stats)[1], 3),
                 dimnames = list(combs, names(stats(rnorm(100))),
                                 c('dist', 'sp', 'dupl')))
    }
    a[i, , 'dist'] <- sdat_stats[, 1]
    a[i, , 'sp'] <- sdat_stats[, 2]
    a[i, , 'dupl'] <- sdat_stats[, 3]
  }
}

stats_df <- data.frame(a[, , 1:3])
stats_df <- na.omit(stats_df)

ggplot(stats_df, aes(Median.dist)) +
  geom_density()

ggplot(stats_df, aes(Mean.dist)) +
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

ggplot(stats_df, aes(Mean.dist, Mean.sp)) +
  geom_point()

ggplot(stats_df, aes(Median.dist, Median.sp)) +
  geom_point()

ggplot(stats_df, aes(Mean.dist, Mean.dupl)) +
  geom_point()

ggplot(stats_df, aes(Median.dist, Median.dupl)) +
  geom_point()
