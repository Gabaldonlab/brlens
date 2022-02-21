# Random phylome stats
# Moisès Bernabeu
# Barcelona, February 2022

# Pakcages loading ----
library(ggpubr)
library(ggplot2)
theme_set(theme_bw())

# Functions definitions ----
get_other <- function(x, ref) {
  if (x['from'] != ref) {
    y <- x['from']
  } else {
    y <- x['to']
  }
  return(y)
}

median_sp <- function(x) {
  by(x, spdat[, 'species_to'],
     median, na.rm = TRUE)
}

# Analysis ----
dat <- read.csv('../outputs/rand_phylome_dist.csv')

spdat <- dat[which(dat$mrca_type == 'S' & dat$from == 'A' | 
                     dat$to == 'A' & dat$from != dat$to), ]

spdat$species_to <- apply(spdat, 1, get_other, ref = 'A')

s2s_med_dist <- aggregate(spdat, by = list(spdat$id), FUN = median)[, c(2, 5)]
row.names(s2s_med_dist) <- s2s_med_dist[, 1]
spdat$dist_norm_s2s <- spdat$dist / s2s_med_dist[spdat$id, 2]

med.df <- data.frame(apply(spdat[, c(4:13, 18)], 2, FUN = median_sp))
med.df <- cbind('species_to' = row.names(med.df), med.df)

# Basic descriptive plots
dist.dens <- ggplot(spdat, aes(dist, col = species_to, fill = species_to)) +
  geom_density(alpha = 0.6) +
  xlim(0, 5) +
  labs(title = 'A to sp.') +
  geom_vline(data = med.df, aes(xintercept = dist,
                                color = species_to),
             show.legend = FALSE)

dist.dens.mrca <- ggplot(spdat, aes(dist_norm_mrca, col = species_to,
                                    fill = species_to)) +
  geom_density(alpha = 0.6) +
  xlim(0, 5) +
  labs(title = 'A to sp. mrca norm')+
  geom_vline(data = med.df, aes(xintercept = dist_norm_mrca,
                                color = species_to),
             show.legend = FALSE)

dist.dens.st <- ggplot(spdat, aes(dist_norm_st, col = species_to,
                                  fill = species_to)) +
  geom_density(alpha = 0.6) +
  xlim(0, 5) +
  labs(title = 'A to sp. subtree norm') +
  geom_vline(data = med.df, aes(xintercept = dist_norm_st,
                                color = species_to),
             show.legend = FALSE)

dist.dens.width <- ggplot(spdat, aes(dist_norm_width, col = species_to,
                                     fill = species_to)) +
  geom_density(alpha = 0.6) +
  xlim(0, 5) +
  labs(title = 'A to sp. tree width norm') +
  geom_vline(data = med.df, aes(xintercept = dist_norm_width,
                                color = species_to),
             show.legend = FALSE)

dist.dens.root <- ggplot(spdat, aes(dist_norm_root, col = species_to,
                                    fill = species_to)) +
  geom_density(alpha = 0.6) +
  xlim(0, 5) +
  labs(title = 'A to sp. root-to-tip norm') +
  geom_vline(data = med.df, aes(xintercept = dist_norm_root,
                                color = species_to),
             show.legend = FALSE)

dist.dens.s2s <- ggplot(spdat, aes(dist_norm_s2s, col = species_to,
                                   fill =  species_to)) +
  geom_density(alpha = 0.6) +
  xlim(0, 5) +
  labs(title = 'A to sp. seq-to-seq norm') +
  geom_vline(data = med.df, aes(xintercept = dist_norm_s2s,
                                color = species_to),
             show.legend = FALSE)

# pdf('../outputs/dist_dens.pdf', width = 10, height = 6)
ggarrange(dist.dens, dist.dens.width, dist.dens.root, dist.dens.st,
          dist.dens.mrca, dist.dens.s2s, align = 'hv', common.legend = TRUE,
          legend = 'bottom')
# dev.off()

# pdf('../outputs/dist_dens_sep.pdf', width = 10, height = 6)
ggplot(spdat, aes(dist, col = species_to, fill = species_to)) +
  geom_density(alpha = 0.6, show.legend = FALSE) +
  geom_vline(data = med.df, aes(xintercept = dist), lty = 4) +
  facet_wrap(~species_to, scales = 'free') +
  xlim(0, 5) +
  labs(title = 'A to sp. raw distance')

ggplot(spdat, aes(dist_norm_width, col = species_to, fill = species_to)) +
  geom_density(alpha = 0.6, show.legend = FALSE) +
  facet_wrap(~species_to, scales = 'free') +
  geom_vline(data = med.df, aes(xintercept = dist_norm_width), lty = 4) +
  xlim(0, 5) +
  labs(title = 'A to sp. tree width normalised distance')

ggplot(spdat, aes(dist_norm_root, col = species_to, fill = species_to)) +
  geom_density(alpha = 0.6, show.legend = FALSE) +
  geom_vline(data = med.df, aes(xintercept = dist_norm_root), lty = 4) +
  facet_wrap(~species_to, scales = 'free') +
  xlim(0, 5) +
  labs(title = 'A to sp. root-to-tip normalized distance')

ggplot(spdat, aes(dist_norm_st, col = species_to, fill = species_to)) +
  geom_density(alpha = 0.6, show.legend = FALSE) +
  geom_vline(data = med.df, aes(xintercept = dist_norm_st), lty = 4) +
  facet_wrap(~species_to, scales = 'free') +
  xlim(0, 5) +
  labs(title = 'A to sp. subtree normalised distance')

ggplot(spdat, aes(dist_norm_mrca, col = species_to, fill = species_to)) +
  geom_density(alpha = 0.6, show.legend = FALSE) +
  geom_vline(data = med.df, aes(xintercept = dist_norm_mrca), lty = 4) +
  facet_wrap(~species_to, scales = 'free') +
  xlim(0, 5) +
  labs(title = 'A to sp. MRCA paris normalised distance')

ggplot(spdat, aes(dist_norm_s2s, col = species_to, fill = species_to)) +
  geom_density(alpha = 0.6, show.legend = FALSE) +
  geom_vline(data = med.df, aes(xintercept = dist_norm_s2s), lty = 4) +
  facet_wrap(~species_to, scales = 'free') +
  xlim(0, 5) +
  labs(title = 'A to sp. seq to seq normalised distance')
# dev.off()

A_sort <- med.df[order(med.df$dist_norm_mrca), 'species_to']

sp_dat <- read.csv('../outputs/sptree_dist.csv')
a <- med.df[, -1] - spdat[1:4, ]$dist
apply(a[-4, ], 2, sum)

write.csv(a, file = '../outputs/norm_tau_differences.csv')
