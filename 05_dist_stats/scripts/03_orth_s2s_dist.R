# Data analysis of phylome distances
# Moisès Bernabeu
# Barcelona, February, 2022

# Libraries loading ----
library(e1071)
library(corrplot)

library(ggplot2)
library(ggpubr)
theme_set(theme_bw())

# Saccharomyces phylome ----
dat <- read.csv('../data/0005_dist_noh.csv')
spdat <- dat[which(dat$mrca_type == 'S' & dat$from_sp == 'YEAST'), ]

# Basic descriptive plots
dist.dens <- ggplot(spdat, aes(dist, col = to_sp, fill = to_sp)) +
  geom_density(alpha = 0.6) +
  xlim(0, 5) +
  labs(title = 'YEAST to sp.')

ndist.dens <- ggplot(spdat, aes(dist_norm, col = to_sp, fill = to_sp)) +
  geom_density(alpha = 0.6) +
  xlim(0, 5) +
  labs(title = 'YEAST to sp.')

# pdf('../outputs/0005_dist_dens.pdf', width = 9, height = 4)
ggarrange(dist.dens, ndist.dens, align = 'h', common.legend = TRUE,
          legend = 'bottom')
# dev.off()

# pdf('../outputs/0005_dist_dens_sep.pdf', width = 10, height = 6)
ggplot(spdat, aes(dist, col = to_sp, fill = to_sp)) +
  geom_density(alpha = 0.6, show.legend = FALSE) +
  facet_wrap(~to_sp, scales = 'free') +
  xlim(0, 5) +
  labs(title = 'Yeast to sp. raw distance')
# dev.off()

# pdf('../outputs/0005_ndist_dens_sep.pdf', width = 10, height = 6)
ggplot(spdat, aes(dist_norm, col = to_sp, fill = to_sp)) +
  geom_density(alpha = 0.6, show.legend = FALSE) +
  facet_wrap(~to_sp, scales = 'free') +
  xlim(0, 5) +
  labs(title = 'Yeast to sp. "normalized" distance')
# dev.off()

# Paired plots
a <- ggplot(spdat, aes(dist, sp)) +
  geom_point()
b <- ggplot(spdat, aes(dist, dupl)) +
  geom_point()
c <- ggplot(spdat, aes(dist, dupl / sp)) +
  geom_point()

# pdf('../outputs/0005_dist_pairs.pdf', width = 12, height = 3)
ggarrange(a, b, c, hjust = 'h', nrow = 1)
# dev.off()

# Extract protein codes for rare peaks in distribution
prots <- spdat$prot[which(spdat$to_sp %in% c('CANAL', 'CANTR', 'CANDU',
                                             'CLALS', 'DEBHA', 'LODEL',
                                             'PICGU', 'PICST') &
                            spdat$dist < 1)]

names(which(table(prots) >= 8))
names(which(table(prots) >= 7))

# Human phylome ----
dat <- read.csv('../data/0076_dist_noh.csv')
spdat <- dat[which(dat$mrca_type == 'S' & dat$from_sp == 'HUMAN'), ]

# Basic descriptive plots
dist.dens <- ggplot(spdat, aes(dist, col = to_sp, fill = to_sp)) +
  geom_density(alpha = 0.6) +
  xlim(0, 1) +
  labs(title = 'HUMAN to sp.')

ndist.dens <- ggplot(spdat, aes(dist_norm, col = to_sp, fill = to_sp)) +
  geom_density(alpha = 0.6) +
  xlim(0, 1) +
  labs(title = 'HUMAN to sp.')

# pdf('../outputs/0076_dist_dens.pdf', width = 9, height = 4)
ggarrange(dist.dens, ndist.dens, align = 'h', common.legend = TRUE,
          legend = 'bottom')
# dev.off()

# pdf('../outputs/0076_dist_dens_sep.pdf', width = 10, height = 6)
ggplot(spdat, aes(dist, col = to_sp, fill = to_sp)) +
  geom_density(alpha = 0.6, show.legend = FALSE) +
  facet_wrap(~to_sp, scales = 'free') +
  # xlim(0, 50) +
  labs(title = 'Human to sp. raw distance')
# dev.off()

# pdf('../outputs/0076_dist_dens_sep.pdf', width = 10, height = 6)
ggplot(spdat, aes(dist_norm, col = to_sp, fill = to_sp)) +
  geom_density(alpha = 0.6, show.legend = FALSE) +
  facet_wrap(~to_sp, scales = 'free') +
  xlim(0, 100) +
  labs(title = 'Human to sp. "normalized" distance')
# dev.off()

# Paired plots
a <- ggplot(spdat, aes(dist, sp)) +
  geom_point()
b <- ggplot(spdat, aes(dist, dupl)) +
  geom_point()
c <- ggplot(spdat, aes(dist, dupl / sp)) +
  geom_point()

# pdf('../outputs/0076_dist_pairs.pdf', width = 12, height = 3)
ggarrange(a, b, c, hjust = 'h', nrow = 1)
# dev.off()
