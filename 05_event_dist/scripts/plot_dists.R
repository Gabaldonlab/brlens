library(ggplot2)
theme_set(theme_bw())

dat <- read.csv('../outputs/0076_dist.csv')

dat <- dat[-which(dat$vert_dist == dat$met_dist), ]

ggplot(dat) +
  geom_density(aes(x = vert_ndist, col = 'vertebrates')) +
  geom_density(aes(x = met_ndist, col = 'metazoan')) +
  xlim(0, quantile(c(dat$met_ndist, dat$vert_ndist), 0.9))
