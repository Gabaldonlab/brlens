# Ploting some tree statistics
# Moisès Bernabeu
# Barcelona, May 2022

library(ggplot2)
library(ggpubr)

theme_set(theme_bw())

dat <- read.csv('../outputs/tree_stats.csv', sep = '\t')

ggplot(dat, aes(tlen_leafno_ratio, mean_brlens, colour = duprate)) +
  geom_point() +
  xlab('Tree length / leaf no.') +
  ylab('Mean branch length')

ggplot(dat, aes(tlen_leafno_ratio, median_brlens, colour = duprate)) +
  geom_point() +
  xlab('Tree length / leaf no.') +
  ylab('Median branch length')

ggplot(dat, aes(tlen_leafno_ratio, median_r2t, colour = duprate)) +
  geom_point() +
  xlab('Tree length / leaf no.') +
  ylab('Median root-to-tip distance')

a <- ggplot(dat, aes(median_r2t, mean_r2t, colour = duprate)) +
  geom_point(alpha = 0.6) +
  xlim(0, 80) +
  ylim(0, 80) +
  geom_smooth(method = 'lm') +
  xlab('Median root-to-tip distances') +
  ylab('Mean root-to-tip distances') +
  scale_colour_continuous(high = "#132B43", low = "#56B1F7")

b <- ggplot(dat, aes(median_brlens, mean_brlens, colour = duprate)) +
  geom_point(alpha = 0.6) +
  xlim(0, 0.75) +
  ylim(0, 2) +
  xlab('Median branch length') +
  ylab('Mean branch length') +
  scale_colour_continuous(high = "#132B43", low = "#56B1F7")

ggarrange(a, b, labels = c('a)', 'b)'), common.legend = TRUE)

dat <- read.csv('../../01_distances/data/trees_stats.csv')
dat <- na.omit(dat)

ggplot(dat, aes(norm_tlen_leafno_ratio)) +
  geom_density()

a <- ggplot(dat, aes(norm_median_r2t, tree_tlen_leafno_ratio, colour = tree_duprate)) +
  geom_point(alpha = 0.6) +
  scale_colour_continuous(high = "#132B43", low = "#56B1F7") +
  xlab('Normalising factor') +
  ylab('Tree length / leaf no.')

b <- ggplot(dat, aes(norm_median_r2t, norm_tlen_leafno_ratio, colour = tree_duprate)) +
  geom_point(alpha = 0.6) +
  scale_colour_continuous(high = "#132B43", low = "#56B1F7") +
  xlab('Normalising factor') +
  ylab('Primates subtree length / leaf no.')

ggarrange(a, b, labels = c('a)', 'b)'), common.legend = TRUE)

ggplot(dat, aes(tree_mean_brlens, tree_tlen_leafno_ratio, colour = tree_duprate)) +
  geom_point() +
  xlab('Mean branch length') +
  ylab('Tree length / leaf no.') +
  labs(colour = 'Dup. rate') +
  scale_colour_continuous(high = "#132B43", low = "#56B1F7")

a <- ggplot(dat, aes(tree_skew_r2t)) +
  geom_density(fill = 'steelblue', colour = 'steelblue', alpha = 0.6) +
  ylab('Density') +
  xlab('Root-to-tip distances skewness')
  

b <- ggplot(dat, aes(tree_skew_brlens)) +
  geom_density(fill = 'steelblue', colour = 'steelblue', alpha = 0.6) +
  ylab('Density') +
  xlab('Branch lengths skewness')

ggarrange(a, b, labels = c('a)', 'b)'), common.legend = TRUE)
