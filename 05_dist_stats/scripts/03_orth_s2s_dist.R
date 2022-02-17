# Data analysis of phylome distances
# Moisès Bernabeu
# Barcelona, February, 2022

# Libraries loading ----
library(e1071)
library(corrplot)

library(ggplot2)
library(ggpubr)
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

median_sp <- function(x) {
  by(x, spdat[, 'species_to'],
     median, na.rm = TRUE)
}

# Saccharomyces phylome ----
dat <- read.csv('../data/0005_dists_noh.csv')
spdat <- dat[which(dat$mrca_type == 'S' & dat$from_sp == 'YEAST' | 
                     dat$to_sp == 'YEAST' & dat$from_sp != dat$to_sp), ]

spdat$species_to <- apply(spdat, 1, get_other, ref = 'YEAST')

# Basic descriptive plots
dist.dens <- ggplot(spdat, aes(dist, col = species_to, fill = species_to)) +
  geom_density(alpha = 0.6) +
  xlim(0, 5) +
  labs(title = 'YEAST to sp.')

dist.dens.mrca <- ggplot(spdat, aes(dist_norm_mrca, col = species_to,
                                    fill = species_to)) +
  geom_density(alpha = 0.6) +
  xlim(0, 5) +
  labs(title = 'YEAST to sp. mrca norm')

dist.dens.st <- ggplot(spdat, aes(dist_norm_st, col = species_to,
                                  fill = species_to)) +
  geom_density(alpha = 0.6) +
  xlim(0, 5) +
  labs(title = 'YEAST to sp. subtree norm')

dist.dens.width <- ggplot(spdat, aes(dist_norm_width, col = species_to,
                                     fill = species_to)) +
  geom_density(alpha = 0.6) +
  xlim(0, 5) +
  labs(title = 'YEAST to sp. tree width norm')

dist.dens.root <- ggplot(spdat, aes(dist_norm_root, col = species_to,
                                    fill = species_to)) +
  geom_density(alpha = 0.6) +
  xlim(0, 5) +
  labs(title = 'YEAST to sp. root-to-tip norm')

# pdf('../outputs/0005_dist_dens.pdf', width = 10, height = 6)
ggarrange(dist.dens, dist.dens.width, dist.dens.root, dist.dens.st,
          dist.dens.mrca, align = 'hv', common.legend = TRUE,
          legend = 'bottom')
# dev.off()

med.df <- data.frame(apply(spdat[, 8:17], 2, FUN = median_sp))
med.df <- cbind('species_to' = row.names(med.df), med.df)

ggplot(spdat, aes(dist_norm_mrca, col = species_to, fill = species_to)) +
  geom_density(alpha = 0.4) +
  xlim(0, 5) +
  labs(title = 'YEAST to sp. tree width norm') +
  geom_vline(data = med.df, aes(xintercept = dist_norm_mrca,
                                color = species_to),
             show.legend = FALSE)

yeast_sort <- med.df[order(med.df$dist_norm_mrca), 'species_to']

# pdf('../outputs/0005_dist_dens_sep.pdf', width = 10, height = 6)
ggplot(spdat, aes(dist, col = species_to, fill = species_to)) +
  geom_density(alpha = 0.6, show.legend = FALSE) +
  geom_vline(data = med.df, aes(xintercept = dist), lty = 4) +
  facet_wrap(~species_to, scales = 'free') +
  xlim(0, 5) +
  labs(title = 'Yeast to sp. raw distance')

ggplot(spdat, aes(dist_norm_width, col = species_to, fill = species_to)) +
  geom_density(alpha = 0.6, show.legend = FALSE) +
  facet_wrap(~species_to, scales = 'free') +
  geom_vline(data = med.df, aes(xintercept = dist_norm_width), lty = 4) +
  xlim(0, 5) +
  labs(title = 'Yeast to sp. tree width normalised distance')

ggplot(spdat, aes(dist_norm_root, col = species_to, fill = species_to)) +
  geom_density(alpha = 0.6, show.legend = FALSE) +
  geom_vline(data = med.df, aes(xintercept = dist_norm_root), lty = 4) +
  facet_wrap(~species_to, scales = 'free') +
  xlim(0, 5) +
  labs(title = 'Yeast to sp. root-to-tip normalized distance')

ggplot(spdat, aes(dist_norm_st, col = species_to, fill = species_to)) +
  geom_density(alpha = 0.6, show.legend = FALSE) +
  geom_vline(data = med.df, aes(xintercept = dist_norm_st), lty = 4) +
  facet_wrap(~species_to, scales = 'free') +
  xlim(0, 5) +
  labs(title = 'Yeast to sp. subtree normalised distance')

ggplot(spdat, aes(dist_norm_mrca, col = species_to, fill = species_to)) +
  geom_density(alpha = 0.6, show.legend = FALSE) +
  geom_vline(data = med.df, aes(xintercept = dist_norm_mrca), lty = 4) +
  facet_wrap(~species_to, scales = 'free') +
  xlim(0, 5) +
  labs(title = 'Yeast to sp. MRCA paris normalised distance')
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
length(prots)
length(table(prots))
length(names(which(table(prots) >= 8)))
names(which(table(prots) >= 7))

# Human phylome ----
dat <- read.csv('../data/0076_dists_noh.csv')
spdat <- dat[which(dat$mrca_type == 'S' & dat$from_sp == 'HUMAN' | 
                     dat$to_sp == 'HUMAN' & dat$from_sp != dat$to_sp), ]

spdat$species_to <- apply(spdat, 1, get_other, ref = 'HUMAN')

# Basic descriptive plots
dist.dens <- ggplot(spdat, aes(dist, col = species_to, fill = species_to)) +
  geom_density(alpha = 0.6) +
  xlim(0, 5) +
  labs(title = 'HUMAN to sp.')

dist.dens.mrca <- ggplot(spdat, aes(dist_norm_mrca, col = species_to,
                                    fill = species_to)) +
  geom_density(alpha = 0.6) +
  xlim(0, 5) +
  labs(title = 'HUMAN to sp. mrca norm')

dist.dens.st <- ggplot(spdat, aes(dist_norm_st, col = species_to,
                                  fill = species_to)) +
  geom_density(alpha = 0.6) +
  xlim(0, 5) +
  labs(title = 'HUMAN to sp. subtree norm')

dist.dens.width <- ggplot(spdat, aes(dist_norm_width, col = species_to,
                                     fill = species_to)) +
  geom_density(alpha = 0.6) +
  xlim(0, 5) +
  labs(title = 'HUMAN to sp. tree width norm')

dist.dens.root <- ggplot(spdat, aes(dist_norm_root, col = species_to,
                                    fill = species_to)) +
  geom_density(alpha = 0.6) +
  xlim(0, 5) +
  labs(title = 'HUMAN to sp. root-to-tip norm')

# pdf('../outputs/0076_dist_dens.pdf', width = 10, height = 6)
ggarrange(dist.dens, dist.dens.width, dist.dens.root, dist.dens.st,
          dist.dens.mrca, align = 'hv', common.legend = TRUE,
          legend = 'bottom')
# dev.off()

med.df <- data.frame(apply(spdat[, 8:17], 2, FUN = median_sp))
med.df <- cbind('species_to' = row.names(med.df), med.df)

ggplot(spdat, aes(dist_norm_mrca, col = species_to, fill = species_to)) +
  geom_density(alpha = 0.4) +
  xlim(0, 5) +
  labs(title = 'HUMAN to sp. tree width norm') +
  geom_vline(data = med.df, aes(xintercept = dist_norm_mrca,
                                color = species_to),
             show.legend = FALSE)

human_sort <- med.df[order(med.df$dist_norm_mrca), 'species_to']

# pdf('../outputs/0076_dist_dens_sep.pdf', width = 10, height = 6)
ggplot(spdat, aes(dist, col = species_to, fill = species_to)) +
  geom_density(alpha = 0.6, show.legend = FALSE) +
  geom_vline(data = med.df, aes(xintercept = dist), lty = 4) +
  facet_wrap(~species_to, scales = 'free') +
  xlim(0, 7) +
  labs(title = 'HUMAN to sp. raw distance')

ggplot(spdat, aes(dist_norm_width, col = species_to, fill = species_to)) +
  geom_density(alpha = 0.6, show.legend = FALSE) +
  facet_wrap(~species_to, scales = 'free') +
  geom_vline(data = med.df, aes(xintercept = dist_norm_width), lty = 4) +
  xlim(0, 7) +
  labs(title = 'HUMAN to sp. tree width normalised distance')

ggplot(spdat, aes(dist_norm_root, col = species_to, fill = species_to)) +
  geom_density(alpha = 0.6, show.legend = FALSE) +
  geom_vline(data = med.df, aes(xintercept = dist_norm_root), lty = 4) +
  facet_wrap(~species_to, scales = 'free') +
  xlim(0, 7) +
  labs(title = 'HUMAN to sp. root-to-tip normalized distance')

ggplot(spdat, aes(dist_norm_st, col = species_to, fill = species_to)) +
  geom_density(alpha = 0.6, show.legend = FALSE) +
  geom_vline(data = med.df, aes(xintercept = dist_norm_st), lty = 4) +
  facet_wrap(~species_to, scales = 'free') +
  xlim(0, 7) +
  labs(title = 'HUMAN to sp. subtree normalised distance')

ggplot(spdat, aes(dist_norm_mrca, col = species_to, fill = species_to)) +
  geom_density(alpha = 0.6, show.legend = FALSE) +
  geom_vline(data = med.df, aes(xintercept = dist_norm_mrca), lty = 4) +
  facet_wrap(~species_to, scales = 'free') +
  xlim(0, 7) +
  labs(title = 'HUMAN to sp. MRCA paris normalised distance')
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

# Species trees information ----
phy5_dat <- read.csv('../data/0005_sptree_sptree_dist.csv')

phy5_spdat <- phy5_dat[which(phy5_dat$from_sp == 'YEAST' |
                                 phy5_dat$to_sp == 'YEAST'), ]
phy5_spdat$species_to <- apply(phy5_spdat, 1, get_other, ref = 'YEAST')

yeast_sort_sptree <- phy5_spdat[order(phy5_spdat$dist_norm_mrca),
                                'species_to']

yeast_sort_comp <- cbind('Phylome' = yeast_sort,
                         'Species tree' = yeast_sort_sptree,
                         'Equal' = yeast_sort == yeast_sort_sptree)

write.csv(yeast_sort_comp, file = '../outputs/0005_sort_comp.csv')

phy76_dat <- read.csv('../data/0076_sptree_sptree_dist.csv')

phy76_spdat <- phy76_dat[which(phy76_dat$from_sp == 'HUMAN' |
                                 phy76_dat$to_sp == 'HUMAN'), ]
phy76_spdat$species_to <- apply(phy76_spdat, 1, get_other, ref = 'HUMAN')

human_sort_sptree <- phy76_spdat[order(phy76_spdat$dist_norm_mrca),
                                 'species_to']

human_sort_comp <- cbind('Phylome' = human_sort,
                         'Species tree' = human_sort_sptree,
                         'Equal' = human_sort == human_sort_sptree)

write.csv(human_sort_comp, file = '../outputs/0076_sort_comp.csv')
