# Plots MSc Thesis
# Moisès Bernabeu
# Barcelona, April 2022

library(treeio)
library(ggtree)

library(ggridges)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(gridExtra)

theme_set(theme_bw())

# Definitions ----
get_other <- function(x, ref) {
  if (x['from_sp'] != ref) {
    y <- x['from_sp']
  } else {
    y <- x['to_sp']
  }
  return(y)
}

# Import and save image of data ----

# Read raw data
ryedat <- read.csv('../../02_get_distances/outputs/0005_dist.csv')
rhudat <- read.csv('../../02_get_distances/outputs/0076_dist.csv')

yespt <- read.csv('../../02_get_distances/outputs/0005_sptree_dist.csv')
huspt <- read.csv('../../02_get_distances/outputs/0076_sptree_dist.csv')

# Get species to
yedat <- ryedat
yedat$sp_to <- apply(ryedat, 1, get_other, ref = 'YEAST')
hudat <- rhudat
hudat$sp_to <- apply(rhudat, 1, get_other, ref = 'HUMAN')

yespt$sp_to <- apply(yespt, 1, get_other, ref = 'YEAST')
huspt$sp_to <- apply(huspt, 1, get_other, ref = 'HUMAN')

# Filter data
yedat <- yedat[which(yedat$mrca_type == 'S' & yedat$sp_to != 'YEAST' &
                        (yedat$to_sp == 'YEAST' | yedat$from_sp == 'YEAST') &
                       yedat$dist != 0), ]
hudat <- hudat[which(hudat$mrca_type == 'S' & hudat$sp_to != 'HUMAN' &
                       (hudat$to_sp == 'HUMAN' | hudat$from_sp == 'HUMAN') &
                       hudat$dist != 0), ]

yespt <- yespt[which(yespt$mrca_type == 'S' & yespt$sp_to != 'YEAST' &
                       (yespt$sp_to == 'YEAST' | yespt$from_sp == 'YEAST')), ]
huspt <- huspt[which(huspt$mrca_type == 'S' & huspt$sp_to != 'HUMAN' &
                       (huspt$sp_to == 'HUMAN' | huspt$from_sp == 'HUMAN')), ]


# save(hudat, yedat, yespt, huspt, file = '../data/seed2sp_dist.Rdata')

# Plots ----
load('../data/seed2sp_dist.Rdata')

# Example trees
set.seed(12031963)

tree <- rtree(20, tip.label = LETTERS[1:20])

d1920 <- c(20, 37, 38, 39, 19)
d0914 <- c(14, 9, 30, 28, 27, 26, 22, 21, 31, 32, 35)
noselnodes <- noselnodes <- c(1:39)[-c(d1920, d0914)]

cols <- c(rep('d(A, N)', length(d1920)),
          rep('d(C, P)', length(d0914)),
          rep(NA, length(noselnodes)))

d <- data.frame(node = c(d1920, d0914, noselnodes), group = cols)

tplot <- ggtree(tree, size = 1, aes(col = group)) %<+% d +
  geom_tippoint(aes(subset = node %in% c(19, 20)), size = 2, col = 'darkorange3') +
  geom_tippoint(aes(subset = node %in% c(9, 14)), size = 2, col = 'steelblue') +
  geom_tiplab(offset = 0.04, show.legend = FALSE) +
  # geom_tiplab(offset = 0.5, aes(label = node)) +
  scale_color_manual(values = c('d(C, P)' = 'steelblue',
                                'd(A, N)' = 'darkorange3'),
                     na.value = 'black') +
  labs(color = 'Distances')

# pdf('../msct_plots/tree.pdf', width = 4, height = 5.5)
tplot
# dev.off()

# Density plots with ggridges
yedens <- yedat %>%
  mutate(sp_to = fct_reorder(.f = sp_to, .x = dist, .fun = median)) %>%
  ggplot(aes(x = dist, y = sp_to, fill = sp_to)) +
  geom_density_ridges2(show.legend = FALSE) +
  xlim(0, 5) +
  xlab('Distance to S. cerevisiae') +
  ylab('Density')

# pdf('../msct_plots/ggridg_seed2sp_yeast.pdf', width = 2.5, height = 7)
yedens
# dev.off()

hudens <- hudat %>%
  mutate(sp_to = fct_reorder(.f = sp_to, .x = dist, .fun = median)) %>%
  ggplot(aes(x = dist, y = sp_to, fill = sp_to)) +
  geom_density_ridges2(show.legend = FALSE) +
  xlim(0, 10) +
  xlab('Distance to H. sapiens') +
  ylab('Density')

# pdf('../msct_plots/ggridg_seed2sp_human.pdf', width = 2.5, height = 7)
hudens
# dev.off()

# pdf('../msct_plots/ggridg_seed2sp_both.pdf', width = 5, height = 7)
ggarrange(yedens, hudens, labels = 'auto', align = 'hv')
# dev.off()

y <- read.csv('../../02_get_distances/outputs/0005_sptree_dist.csv')
y <- y[, c(4, 6, 10)]

nsp <- sum(!duplicated(y$from_sp))
sps <- y$from_sp[!duplicated(y$from_sp)]

m <- matrix(NA, nrow = nsp, ncol = nsp)
for (i in 1:nsp) {
  for (j in i:nsp) {
    if (sps[i] != sps[j]) {
      m[i, j] <- y[which((y$from_sp == sps[i] & y$to_sp == sps[j]) |
                           (y$from_sp == sps[j] & y$to_sp == sps[i])), 3]
      m[j, i] <- y[which((y$from_sp == sps[i] & y$to_sp == sps[j]) |
                           (y$from_sp == sps[j] & y$to_sp == sps[i])), 3]
    }
  }
}

# pdf('../msct_plots/heatmap.pdf', width = 5, height = 5)
heatmap(m, Rowv = NA, Colv = NA,
        labRow = LETTERS[1:20],
        labCol = LETTERS[1:20])
# dev.off()

tplots <- c()
for (i in 1:20) {
  n <- rpois(1, 20)
  tree <- rtree(n)
  tplots[[i]] <- ggtree(tree)
}

# pdf('../msct_plots/phylome_2.pdf', width = 4, height = 7)
do.call(grid.arrange, tplots)
# dev.off()

# Density plots with ggridges of normalised distances
yendens <- yedat %>%
  mutate(sp_to = fct_reorder(.f = sp_to, .x = ndist_A, .fun = median)) %>%
  ggplot(aes(x = ndist_A, y = sp_to, fill = sp_to)) +
  geom_density_ridges2(show.legend = FALSE) +
  xlim(0, 5) +
  xlab('Distance to S. cerevisiae') +
  ylab('Density')

# pdf('../msct_plots/ggridg_seed2sp_n_yeast.pdf', width = 2.5, height = 7)
yendens
# dev.off()

hundens <- hudat %>%
  mutate(sp_to = fct_reorder(.f = sp_to, .x = ndist_A, .fun = median)) %>%
  ggplot(aes(x = ndist_A, y = sp_to, fill = sp_to)) +
  geom_density_ridges2(show.legend = FALSE) +
  xlim(0, 10) +
  xlab('Distance to H. sapiens') +
  ylab('Density')

# pdf('../msct_plots/ggridg_seed2sp_n_human.pdf', width = 2.5, height = 7)
hundens
# dev.off()

# pdf('../msct_plots/ggridg_seed2sp_n_both.pdf', width = 5, height = 7)
ggarrange(yendens, hundens, labels = 'auto', align = 'hv')
# dev.off()

yedplot <- ggplot(yedat, aes(dist, colour = sp_to)) +
  geom_density(show.legend = FALSE) +
  xlim(0, quantile(yedat$dist, 0.999)) +
  ylab('Density') +
  xlab('Raw distance')

yendplot <- ggplot(yedat, aes(ndist_A, colour = sp_to)) +
  geom_density(show.legend = FALSE) +
  xlim(0, quantile(yedat$dist, 0.999)) +
  ylab('Density') +
  xlab('Normalised distance')

# pdf('../outputs/yeast_densities.pdf', width = 6.3, height = 2.3)
ggarrange(yedplot, yendplot, labels = 'auto', align = 'hv')
# dev.off()

hudplot <- ggplot(hudat, aes(dist, colour = sp_to)) +
  geom_density(show.legend = FALSE) +
  xlim(0, quantile(hudat$dist, 0.999)) +
  ylab('Density') +
  xlab('Raw distance')

hundplot <- ggplot(hudat, aes(ndist_A, colour = sp_to)) +
  geom_density(show.legend = TRUE) +
  xlim(0, quantile(hudat$dist, 0.999)) +
  ylab('Density') +
  xlab('Normalised distance')

# pdf('../outputs/human_densities.pdf', width = 6.3, height = 2.3)
ggarrange(hudplot, hundplot, labels = 'auto', align = 'hv')
# dev.off()

