library(tidyr)
library(ggplot2)
theme_set(theme_bw())

dat <- read.csv('../outputs/0076_dist.csv')

ggplot(dat, aes(norm_median, norm_D / norm_S)) +
  geom_point()

gdupl <- gather(dat[, c('max_median', 'min_median')])

ggplot(gdupl, aes(y = value, x = key, col = key, fill = key)) +
  geom_boxplot(alpha = 0.6) +
  ylim(0, quantile(gdupl$value, 0.9))

gdupl <- gather(dat[, c('max_mean', 'min_mean')])

ggplot(gdupl, aes(y = value, x = key, col = key, fill = key)) +
  geom_boxplot(alpha = 0.6) +
  ylim(0, quantile(gdupl$value, 0.9))
