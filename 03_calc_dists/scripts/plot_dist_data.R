library(ggplot2)
library(ggpubr)
theme_set(theme_bw())

setwd('~/Documents/github/brlens/03_calc_dists/outputs/')

dists <- read.csv('0469_dist.tsv', sep = '\t')
str(dists)

hist(dists$seed_dist)
hist(dists$seed_ndist, xlim = c(0, 40), breaks = 100)
hist(dists$og_dist)
hist(dists$og_ndist, xlim = c(0, 40), breaks = 100)

a <- ggplot(dists, aes(x = og_ndist)) +
  geom_density() +
  xlim(0, 20)

b <- ggplot(dists, aes(x = seed_ndist)) +
  geom_density() +
  xlim(0, 20)

c <- ggplot(dists, aes(x = og_dist)) +
  geom_density() +
  xlim(0, 20)

d <- ggplot(dists, aes(x = seed_dist)) +
  geom_density() +
  xlim(0, 20)

ggarrange(a, b, c, d, ncol = 2, nrow = 2, align = 'hv')
