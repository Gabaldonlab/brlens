# Data analysis of phylome distances - Combined plots
# Moisès Bernabeu
# Barcelona, February, 2022

# Libraries loading ----
library(ggplot2)
library(ggpubr)
library(tidyr)
library(ggrepel)

library(ggtree)
library(treeio)

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

b +
  geom_label_repel(aes(label = sp_to), min.segment.length = 0) +
  xlim(1.5, 2.5) +
  ylim(2, 2.5)

plot(hu_sp_vs_phy$dist_sp, hu_sp_vs_phy$brl_ndist, type = 'n')
text(hu_sp_vs_phy$dist_sp, hu_sp_vs_phy$brl_ndist, labels = hu_sp_vs_phy$sp_to)

hu.t <- read.newick('../../07_normalization/data/0076_sptree.nwk',
                    node.label = 'support')
hu.tp <- ggtree(hu.t) %<+% hu_sp_vs_phy +
  geom_tiplab(aes(fill = group), geom = 'label', size = 2) +
  # geom_nodelab(aes(label = support, x = branch)) +
  geom_nodepoint(aes(colour = cut(support, 3), x = branch)) +
  scale_color_grey(na.value = NA, start = 0.8, end = 0) +
  xlim(c(0, 7.6))
hu.tp

ye.t <- read.newick('../../07_normalization/data/0005_sptree.nwk',
                    node.label = 'support')
ye.tp <- ggtree(ye.t) %<+% ye_sp_vs_phy +
  geom_tiplab(aes(fill = group), geom = 'label', size = 3) +
  # geom_nodelab(aes(label = support, x = branch)) +
  geom_nodepoint(aes(colour = cut(support, 3), x = branch)) +
  scale_color_grey(na.value = NA, start = 0.8, end = 0) +
  xlim(c(0, 2))
ye.tp

ye.tp + xlim(0, 3) + hu.tp + xlim(0, 8)
