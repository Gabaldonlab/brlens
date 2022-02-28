# Data analysis of phylome distances - Combined plots
# Moisès Bernabeu
# Barcelona, February, 2022

# Libraries loading ----
library(ggplot2)
library(ggpubr)
library(tidyr)

theme_set(theme_bw())

# Loading data ----
ye_distdf <- read.csv(paste0('../outputs/0005_distdf.csv'),
                   row.names = 1)
ye_sp_vs_phy <- read.csv('../outputs/0005_sp_vs_phylome.csv',
                      row.names = 1)
ye_sp_vs_phy$group <- factor(ye_sp_vs_phy$group)

hu_distdf <- read.csv(paste0('../outputs/0076_distdf.csv'),
                      row.names = 1)
hu_sp_vs_phy <- read.csv('../outputs/0076_sp_vs_phylome.csv',
                         row.names = 1)
hu_sp_vs_phy$group <- factor(hu_sp_vs_phy$group)

# Plots ----
# Joint densities
a <- ggplot(ye_distdf, aes(x = raw_dist, colour = sp_to)) +
  geom_density(show.legend = FALSE) +
  xlim(0, 5) +
  labs(title = 'YEAST to species')

b <- ggplot(ye_distdf, aes(x = brl_ndist, colour = sp_to)) +
  geom_density(show.legend = FALSE) +
  xlim(0, 5) +
  labs(title = 'YEAST to species')

c <- ggplot(hu_distdf, aes(x = raw_dist, colour = sp_to)) +
  geom_density(show.legend = FALSE) +
  xlim(0, 10) +
  labs(title = 'HUMAN to species')

d <- ggplot(hu_distdf, aes(x = brl_ndist, colour = sp_to)) +
  geom_density(show.legend = FALSE) +
  xlim(0, 7) +
  labs(title = 'HUMAN to species')

ggarrange(a, b, c, d, align = 'hv')

a <- ggplot(ye_sp_vs_phy, aes(dist_sp, brl_ndist, colour = group)) +
  geom_point(show.legend = FALSE)

b <- ggplot(hu_sp_vs_phy, aes(dist_sp, brl_ndist, colour = group)) +
  geom_point(show.legend = FALSE)

c <- ggplot(ye_sp_vs_phy, aes(dist_sp, mrca_ndist, colour = group)) +
  geom_point(show.legend = FALSE)

d <- ggplot(hu_sp_vs_phy, aes(dist_sp, mrca_ndist, colour = group)) +
  geom_point(show.legend = FALSE)

ggarrange(a, b, c, d, align = 'hv', labels = 'auto')
