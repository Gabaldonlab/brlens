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
dat <- read.csv('../data/0005_dists_noh.csv')
spdat <- dat[which(dat$mrca_type == 'S' & dat$from_sp == 'YEAST'), ]

# Basic descriptive plots
dist.dens <- ggplot(spdat, aes(dist, col = to_sp, fill = to_sp)) +
  geom_density(alpha = 0.6) +
  xlim(0, 5) +
  labs(title = 'YEAST to sp.')

dist.dens.mrca <- ggplot(spdat, aes(dist_norm_mrca, col = to_sp, fill = to_sp)) +
  geom_density(alpha = 0.6) +
  xlim(0, 5) +
  labs(title = 'YEAST to sp. mrca norm')

dist.dens.st <- ggplot(spdat, aes(dist_norm_st, col = to_sp, fill = to_sp)) +
  geom_density(alpha = 0.6) +
  xlim(0, 5) +
  labs(title = 'YEAST to sp. subtree norm')

dist.dens.width <- ggplot(spdat, aes(dist_norm_width, col = to_sp, fill = to_sp)) +
  geom_density(alpha = 0.6) +
  xlim(0, 5) +
  labs(title = 'YEAST to sp. tree width norm')

dist.dens.root <- ggplot(spdat, aes(dist_norm_root, col = to_sp, fill = to_sp)) +
  geom_density(alpha = 0.6) +
  xlim(0, 5) +
  labs(title = 'YEAST to sp. root-to-tip norm')

# pdf('../outputs/0005_dist_dens.pdf', width = 10, height = 6)
ggarrange(dist.dens, dist.dens.width, dist.dens.root, dist.dens.st,
          dist.dens.mrca, align = 'hv', common.legend = TRUE,
          legend = 'bottom')
# dev.off()

med.df <- data.frame(apply(spdat[, 8:17], 2, FUN = function(x) {by(x, spdat[, 'to_sp'],
                                                        median, na.rm = TRUE)}))
med.df <- cbind('to_sp' = row.names(med.df), med.df)

ggplot(spdat, aes(dist_norm_mrca, col = to_sp, fill = to_sp)) +
  geom_density(alpha = 0.4) +
  xlim(0, 5) +
  labs(title = 'YEAST to sp. tree width norm') +
  geom_vline(data = med.df, aes(xintercept = dist_norm_mrca, color = to_sp),
             show.legend = FALSE)

med.df[order(med.df$dist_norm_mrca), 'to_sp']

# pdf('../outputs/0005_dist_dens_sep.pdf', width = 10, height = 6)
ggplot(spdat, aes(dist, col = to_sp, fill = to_sp)) +
  geom_density(alpha = 0.6, show.legend = FALSE) +
  geom_vline(data = med.df, aes(xintercept = dist), lty = 4) +
  facet_wrap(~to_sp, scales = 'free') +
  xlim(0, 5) +
  labs(title = 'Yeast to sp. raw distance')

ggplot(spdat, aes(dist_norm_width, col = to_sp, fill = to_sp)) +
  geom_density(alpha = 0.6, show.legend = FALSE) +
  facet_wrap(~to_sp, scales = 'free') +
  geom_vline(data = med.df, aes(xintercept = dist_norm_width), lty = 4) +
  xlim(0, 5) +
  labs(title = 'Yeast to sp. tree width normalised distance')

ggplot(spdat, aes(dist_norm_root, col = to_sp, fill = to_sp)) +
  geom_density(alpha = 0.6, show.legend = FALSE) +
  geom_vline(data = med.df, aes(xintercept = dist_norm_root), lty = 4) +
  facet_wrap(~to_sp, scales = 'free') +
  xlim(0, 5) +
  labs(title = 'Yeast to sp. root-to-tip normalized distance')

ggplot(spdat, aes(dist_norm_st, col = to_sp, fill = to_sp)) +
  geom_density(alpha = 0.6, show.legend = FALSE) +
  geom_vline(data = med.df, aes(xintercept = dist_norm_st), lty = 4) +
  facet_wrap(~to_sp, scales = 'free') +
  xlim(0, 5) +
  labs(title = 'Yeast to sp. subtree normalised distance')

ggplot(spdat, aes(dist_norm_mrca, col = to_sp, fill = to_sp)) +
  geom_density(alpha = 0.6, show.legend = FALSE) +
  geom_vline(data = med.df, aes(xintercept = dist_norm_mrca), lty = 4) +
  facet_wrap(~to_sp, scales = 'free') +
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
spdat <- dat[which(dat$mrca_type == 'S' & dat$from_sp == 'HUMAN'), ]

# Basic descriptive plots
dist.dens <- ggplot(spdat, aes(dist, col = to_sp, fill = to_sp)) +
  geom_density(alpha = 0.6) +
  xlim(0, 5) +
  labs(title = 'HUMAN to sp.')

dist.dens.mrca <- ggplot(spdat, aes(dist_norm_mrca, col = to_sp, fill = to_sp)) +
  geom_density(alpha = 0.6) +
  xlim(0, 5) +
  labs(title = 'HUMAN to sp. mrca norm')

dist.dens.st <- ggplot(spdat, aes(dist_norm_st, col = to_sp, fill = to_sp)) +
  geom_density(alpha = 0.6) +
  xlim(0, 5) +
  labs(title = 'HUMAN to sp. subtree norm')

dist.dens.width <- ggplot(spdat, aes(dist_norm_width, col = to_sp, fill = to_sp)) +
  geom_density(alpha = 0.6) +
  xlim(0, 5) +
  labs(title = 'HUMAN to sp. tree width norm')

dist.dens.root <- ggplot(spdat, aes(dist_norm_root, col = to_sp, fill = to_sp)) +
  geom_density(alpha = 0.6) +
  xlim(0, 5) +
  labs(title = 'HUMAN to sp. root-to-tip norm')

# pdf('../outputs/0076_dist_dens.pdf', width = 10, height = 6)
ggarrange(dist.dens, dist.dens.width, dist.dens.root, dist.dens.st,
          dist.dens.mrca, align = 'hv', common.legend = TRUE,
          legend = 'bottom')
# dev.off()

med.df <- data.frame(apply(spdat[, 8:17], 2, FUN = function(x) {by(x, spdat[, 'to_sp'],
                                                                   median, na.rm = TRUE)}))
med.df <- cbind('to_sp' = row.names(med.df), med.df)

ggplot(spdat, aes(dist_norm_mrca, col = to_sp, fill = to_sp)) +
  geom_density(alpha = 0.4) +
  xlim(0, 5) +
  labs(title = 'HUMAN to sp. tree width norm') +
  geom_vline(data = med.df, aes(xintercept = dist_norm_mrca, color = to_sp),
             show.legend = FALSE)

med.df[order(med.df$dist_norm_mrca), 'to_sp']

# pdf('../outputs/0076_dist_dens_sep.pdf', width = 12, height = 8)
ggplot(spdat, aes(dist, col = to_sp, fill = to_sp)) +
  geom_density(alpha = 0.6, show.legend = FALSE) +
  geom_vline(data = med.df, aes(xintercept = dist), lty = 4) +
  facet_wrap(~to_sp, scales = 'free') +
  xlim(0, 5) +
  labs(title = 'HUMAN to sp. raw distance')

ggplot(spdat, aes(dist_norm_width, col = to_sp, fill = to_sp)) +
  geom_density(alpha = 0.6, show.legend = FALSE) +
  facet_wrap(~to_sp, scales = 'free') +
  geom_vline(data = med.df, aes(xintercept = dist_norm_width), lty = 4) +
  xlim(0, 5) +
  labs(title = 'HUMAN to sp. tree width normalised distance')

ggplot(spdat, aes(dist_norm_root, col = to_sp, fill = to_sp)) +
  geom_density(alpha = 0.6, show.legend = FALSE) +
  geom_vline(data = med.df, aes(xintercept = dist_norm_root), lty = 4) +
  facet_wrap(~to_sp, scales = 'free') +
  xlim(0, 5) +
  labs(title = 'HUMAN to sp. root-to-tip normalized distance')

ggplot(spdat, aes(dist_norm_st, col = to_sp, fill = to_sp)) +
  geom_density(alpha = 0.6, show.legend = FALSE) +
  geom_vline(data = med.df, aes(xintercept = dist_norm_st), lty = 4) +
  facet_wrap(~to_sp, scales = 'free') +
  xlim(0, 5) +
  labs(title = 'HUMAN to sp. subtree normalised distance')

ggplot(spdat, aes(dist_norm_mrca, col = to_sp, fill = to_sp)) +
  geom_density(alpha = 0.6, show.legend = FALSE) +
  geom_vline(data = med.df, aes(xintercept = dist_norm_mrca), lty = 4) +
  facet_wrap(~to_sp, scales = 'free') +
  xlim(0, 5) +
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