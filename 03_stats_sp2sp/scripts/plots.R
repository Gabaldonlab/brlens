# Plots for species to species distances
# Moisès Bernabeu
# Barcelona, April 2022

library(ggridges)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(dplyr)

theme_set(theme_bw())

# Plots ----
load('../data/seed2sp_dist.Rdata')

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

yedensind <- yedat %>%
  mutate(sp_to = fct_reorder(.f = sp_to, .x = dist, .fun = median)) %>%
  ggplot(aes(x = dist, fill = sp_to, colour = sp_to)) +
  geom_density(show.legend = FALSE, alpha = 0.6) +
  xlim(0, quantile(yedat$dist, 0.999)) +
  xlab('Distance to S. cerevisiae') +
  ylab('Density') +
  facet_wrap(~sp_to, scales = 'free_y')

yedensind

# Descriptive yeast
# pdf('../msct_plots/yeast_raw_descr.pdf', width = 11.5, height = 5.75)
ggarrange(yedens, yedensind, align = 'v', widths = c(1, 4), labels = 'auto')
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

hudensind <- hudat %>%
  mutate(sp_to = fct_reorder(.f = sp_to, .x = dist, .fun = median)) %>%
  ggplot(aes(x = dist, fill = sp_to, colour = sp_to)) +
  geom_density(show.legend = FALSE, alpha = 0.6) +
  xlim(0, quantile(hudat$dist, 0.999)) +
  xlab('Distance to H. sapiens') +
  ylab('Density') +
  facet_wrap(~sp_to, scales = 'free_y')

hudensind

# Descriptive human
# pdf('../msct_plots/human_raw_descr.pdf', width = 11.5, height = 5.75)
ggarrange(hudens, hudensind, align = 'v', widths = c(1, 4), labels = 'auto')
# dev.off()

# pdf('../msct_plots/ggridg_seed2sp_both.pdf', width = 5, height = 7)
ggarrange(yedens, hudens, labels = 'auto', align = 'hv')
# dev.off()

# Density plots with ggridges of normalised distances
yendens <- yedat %>%
  mutate(sp_to = fct_reorder(.f = sp_to, .x = ndist_A, .fun = median)) %>%
  ggplot(aes(x = ndist_A, y = sp_to, fill = sp_to)) +
  geom_density_ridges2(show.legend = FALSE) +
  xlim(0, 5) +
  xlab('Norm. dist. to S. cerevisiae') +
  ylab('Density')

yendensind <- yedat %>%
  mutate(sp_to = fct_reorder(.f = sp_to, .x = ndist_A, .fun = median)) %>%
  ggplot(aes(x = ndist_A, fill = sp_to, colour = sp_to)) +
  geom_density(show.legend = FALSE, alpha = 0.6) +
  xlim(0, quantile(yedat$ndist_A, 0.99)) +
  xlab('Normalised distance to S. cerevisiae') +
  ylab('Density') +
  facet_wrap(~sp_to, scales = 'free_y')
yendensind

# pdf('../msct_plots/ggridg_seed2sp_n_yeast.pdf', width = 2.5, height = 7)
yendens
# dev.off()

# pdf('../msct_plots/yeast_norm_descr.pdf', width = 11.5, height = 5.75)
ggarrange(yendens, yendensind, align = 'v', widths = c(1, 4), labels = 'auto')
# dev.off()

hundens <- hudat %>%
  mutate(sp_to = fct_reorder(.f = sp_to, .x = ndist_A, .fun = median)) %>%
  ggplot(aes(x = ndist_A, y = sp_to, fill = sp_to)) +
  geom_density_ridges2(show.legend = FALSE) +
  xlim(0, 15) +
  xlab('Norm. dist. to H. sapiens') +
  ylab('Density')

# pdf('../msct_plots/ggridg_seed2sp_n_human.pdf', width = 2.5, height = 7)
hundens
# dev.off()


hundensind <- hudat %>%
  mutate(sp_to = fct_reorder(.f = sp_to, .x = ndist_A, .fun = median)) %>%
  ggplot(aes(x = ndist_A, fill = sp_to, colour = sp_to)) +
  geom_density(show.legend = FALSE, alpha = 0.6) +
  xlim(0, quantile(hudat$ndist_A)) +
  xlab('Normalised distance to H. sapiens') +
  ylab('Density') +
  facet_wrap(~sp_to, scales = 'free_y')

hundensind

# pdf('../msct_plots/human_norm_descr.pdf', width = 11.5, height = 5.75)
ggarrange(hundens, hundensind, align = 'v', widths = c(1, 4), labels = 'auto')
# dev.off()

# pdf('../msct_plots/ggridg_seed2sp_n_both.pdf', width = 5, height = 7)
ggarrange(yendens, hundens, labels = 'auto', align = 'hv')
# dev.off()

# Joint densities ----
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
  ylim(0, 2.5) +
  ylab('Density') +
  xlab('Raw distance')

hundplot <- ggplot(hudat, aes(ndist_A, colour = sp_to)) +
  geom_density(show.legend = FALSE) +
  xlim(0, quantile(hudat$dist, 0.999)) +
  ylab('Density') +
  xlab('Normalised distance') +
  ylim(0, 2.5)

# pdf('../outputs/human_densities.pdf', width = 6.3, height = 2.3)
ggarrange(hudplot, hundplot, labels = 'auto', align = 'hv')
# dev.off()

# pdf('~/Desktop/densplots.pdf', width = 8, height = 5.5)
ggarrange(yedplot, yendplot, hudplot, hundplot, labels = 'auto', align = 'hv')
# dev.off()

yejoined <- ggplot(yedat, aes(ndist_A)) +
  geom_density() +
  xlim(0, quantile(yedat$ndist_A, 0.99)) +
  xlab('Normalised distance') +
  ylab('Density')

ggarrange(yejoined, yendplot)

yedatgr <- yedat %>%
  group_by(sp_to) %>%
  summarise('median_dist' = median(dist), 'mean_dist' = mean(dist),
            'median_ndist' = median(ndist_A), 'mean_ndist' = mean(ndist_A))

yespphy <- merge(yedatgr, yespt[, c('dist', 'ndist_A', 'sp_to')], by = 'sp_to')

hudatgr <- hudat %>%
  group_by(sp_to) %>%
  summarise('median_dist' = median(dist), 'mean_dist' = mean(dist),
            'median_ndist' = median(ndist_A), 'mean_ndist' = mean(ndist_A))

huspphy <- merge(hudatgr, huspt[, c('dist', 'ndist_A', 'sp_to')], by = 'sp_to')

yespphyp <- ggplot(yespphy, aes(dist, median_dist)) +
  geom_point() +
  xlab('Species tree distance') +
  ylab('Phylome median distance') +
  labs(title = 'Yeast phylome')

huspphyp <- ggplot(huspphy, aes(dist, median_dist)) +
  geom_point() +
  xlab('Species tree distance') +
  ylab('Phylome median distance') +
  labs(title = 'Human phylome')

yespphynp <- ggplot(yespphy, aes(dist, median_ndist)) +
  geom_point() +
  xlab('Species tree distance') +
  ylab('Phylome median norm. distance') +
  labs(title = 'Yeast phylome')

huspphynp <- ggplot(huspphy, aes(dist, median_ndist)) +
  geom_point() +
  xlab('Species tree distance') +
  ylab('Phylome median norm. distance') +
  labs(title = 'Human phylome')

# pdf('../msct_plots/spphy_both.pdf', width = 8.4, height = 3.15)
ggarrange(yespphyp, huspphyp, labels = 'auto', align = 'hv')
# dev.off()

# pdf('../msct_plots/spphy_n_both.pdf', width = 8.4, height = 3.15)
ggarrange(yespphynp, huspphynp, labels = 'auto', align = 'hv')
# dev.off()

# pdf('../msct_plots/sptree.pdf', width = 7.5, height = 5.5)
ggarrange(yespphyp, huspphyp, yespphynp, huspphynp,
          labels = 'auto', align = 'hv')
# dev.off()

# pdf('../outputs/normalised_densities.pdf', width = 6.3, height = 2.3)
ggarrange(yendplot + xlab('Yeast to sp. norm. dist.'),
          hundplot + xlab('Human to sp. norm. dist.'),
          align = 'hv', labels = 'auto')
# dev.off()

yedat$from_to <- paste(yedat$from_sp, yedat$to_sp, sep = '-')
yedat_sum <- yedat %>%
  group_by(sp_to) %>%
  summarise('dist' = median(dist),
            'ndist_A' = median(ndist_A, na.rm = TRUE),
            'sp' = mean(sp),
            'dupl' = mean(dupl))

a <- ggplot(yedat_sum, aes(dupl / (sp + dupl), ndist_A)) +
  geom_point() +
  xlab('Duplication rate') +
  ylab('S. cerevisiae to sp. norm. dist.')

hudat$from_to <- paste(hudat$from_sp, hudat$to_sp, sep = '-')
hudat_sum <- hudat %>%
  group_by(sp_to) %>%
  summarise('dist' = median(dist),
            'ndist_A' = median(ndist_A, na.rm = TRUE),
            'sp' = mean(sp),
            'dupl' = mean(dupl))

b <- ggplot(hudat_sum, aes(dupl / (sp + dupl), ndist_A)) +
  geom_point() +
  xlab('Duplication rate') +
  ylab('H. sapiens to sp. norm. dist.')

# pdf('../msct_plots/duprate.pdf', width = 9, height = 3)
ggarrange(a, b, labels = 'auto', align = 'hv')
# dev.off()
