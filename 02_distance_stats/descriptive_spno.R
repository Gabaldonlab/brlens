# Descriptive statistics of species number
# Moisès Bernabeu
# Sant Feliu de Guíxols, September 2022

library(tidyr)
library(dplyr)
library(ggplot2)

tdat <- read.csv('../data/sps_prop.csv')
ddat <- read.csv('../data/mammal_dist.csv')
taxdat <- read.csv('../data/taxa.csv', row.names = 1)

trees <- ddat[which(ddat$tree_leafno - ddat$norm_leafno < 4), 1]
write.csv(trees, row.names = FALSE, file = '../outputs/trees.txt')

primprop <- ddat$norm_leafno / ddat$tree_leafno
names(tdat)

ggplot() +
  geom_density(aes(x = primprop))

props <- gather(tdat, key = 'species', value = 'proportion', -seed)
props$Order <- taxdat[props$species, 'Order']

ordpropr <- props[, -2] %>%
  group_by(seed, Order) %>%
  summarise(sum_prop = sum(proportion, na.rm = TRUE))

ggplot(ordpropr, aes(sum_prop)) +
  geom_density() +
  xlim(0, 1) +
  facet_wrap(~Order, scales = 'free_y') +
  labs(title = 'Tree data') +
  xlab('Proportion of leaves per tree')


