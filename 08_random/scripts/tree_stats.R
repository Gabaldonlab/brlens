# Random phylome stats
# Moisès Bernabeu
# Barcelona, February 2022

# Pakcages loading ----
library(ggpubr)
library(ggplot2)
library(tidyr)
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

spdat <- dat[which(dat$from == 'A' | 
                     dat$to == 'A' & dat$from != dat$to), ]

spdat$species_to <- apply(spdat, 1, get_other, ref = 'A')

s2s_ref <- aggregate(spdat, by = list(spdat$id), FUN = median)[, c(2, 5)]
row.names(s2s_ref) <- s2s_ref[, 1]

mrca_ref <- read.csv('../outputs/mrca.csv', row.names = 1)
rbls_ref <- read.csv('../outputs/rbls.csv', row.names = 1)
root_ref <- read.csv('../outputs/root.csv', row.names = 1)
st_ref <- read.csv('../outputs/st.csv', row.names = 1)

spdat$mrca_dist <- spdat$dist / mrca_ref[spdat$id, 'med']
spdat$root_dist <- spdat$dist / root_ref[spdat$id, 'med']
spdat$st_dist <- spdat$dist / st_ref[spdat$id, 'med']
spdat$s2s_dist <- spdat$dist / s2s_ref[spdat$id, 'dist']

norm_fact1 <- rbls_ref[spdat$id, 'sum_brl'] /
  median(rbls_ref$sum_brl) + rbls_ref[spdat$id, 'twdth'] /
  median(rbls_ref[, 'twdth'])
norm_fact2 <- rbls_ref[spdat$id, 'sum_brl'] /
  (median(rbls_ref$sum_brl))
norm_fact3 <- rbls_ref[spdat$id, 'med_brl'] / median(rbls_ref$med_brl)

x <- 0:400/100
rp <- ggplot() +
  geom_line(aes(x, dgamma(x, 23, 12) + dgamma(x, 1, 3))) +
  ylab('density') +
  xlab('rate')

np <- ggplot() +
  geom_density(aes(x = norm_fact3)) +
  xlab('Median of branches length normalisation factor')

jpeg('')
ggarrange(rp, np, align = 'hv')

spdat$brls1 <- spdat$dist / norm_fact1
spdat$brls2 <- spdat$dist / norm_fact2
spdat$brls3 <- spdat$dist / norm_fact3

x <- 0:400/100
a <- ggplot() +
  geom_line(aes(x, dgamma(x, 23, 12) + dgamma(x, 1, 3))) +
  ylab('density') +
  xlab('rate')
b <- ggplot() +
  geom_density(aes(x = norm_fact3)) +
  xlab('Median of branches length normalisation factor')

rate_norm <- ggarrange(a, b, align = 'h')

spdatg <- gather(spdat[, -c(2, 3)], value = 'value',
                 key = 'key', -c(id, species_to))

apply(spdat[-c(1:3, 5)], 2, quantile)

ggplot(spdatg, aes(x = value, color = species_to)) +
  geom_density() +
  facet_wrap(~key, scales = 'free_y', shrink = TRUE) +
  xlim(0, 10)

a <- ggplot(spdat, aes(x = dist, col = species_to)) +
  geom_density()
b <- ggplot(spdat, aes(x = root_dist, col = species_to)) +
  geom_density() +
  xlim(0, 2.5) +
  xlab('ndist with root to tip median')
c <- ggplot(spdat, aes(x = mrca_dist, col = species_to)) +
  geom_density() +
  xlim(0, 4) +
  xlab('ndist with mrca to tip pairs median')
d <- ggplot(spdat, aes(x = brls3, col = species_to)) +
  geom_density() +
  xlim(0, 8) +
  xlab('ndist with brlens median and width')

ggarrange(a, b, c, d, common.legend = TRUE, align = 'hv', legend = 'bottom')

a <- ggplot(spdat, aes(x = dist, col = species_to)) +
  geom_boxplot()
b <- ggplot(spdat, aes(x = root_dist, col = species_to)) +
  geom_boxplot() +
  xlim(0, 2.5) +
  xlab('ndist with root to tip median')
c <- ggplot(spdat, aes(x = mrca_dist, col = species_to)) +
  geom_boxplot() +
  xlim(0, 4) +
  xlab('ndist with mrca to tip pairs median')
d <- ggplot(spdat, aes(x = brls1, col = species_to)) +
  geom_boxplot() +
  xlim(0, 10) +
  xlab('ndist with standardized brlens sum and width')

ggarrange(a, b, c, d, common.legend = TRUE, align = 'hv')

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

dist.dens.rbls <- ggplot(spdat, aes(dist_norm_rbls, col = species_to,
                                    fill =  species_to)) +
  geom_density(alpha = 0.6) +
  xlim(0, 5) +
  labs(title = 'A to sp. seq-to-seq norm') +
  geom_vline(data = med.df, aes(xintercept = dist_norm_rbls,
                                color = species_to),
             show.legend = FALSE)

# pdf('../outputs/dist_dens.pdf', width = 10, height = 6)
ggarrange(dist.dens, dist.dens.width, dist.dens.root, dist.dens.st,
          dist.dens.mrca, dist.dens.s2s, dist.dens.rbls,
          align = 'hv', common.legend = TRUE,
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

ggplot(spdat, aes(dist_norm_rbls, col = species_to, fill = species_to)) +
  geom_density(alpha = 0.6, show.legend = FALSE) +
  geom_vline(data = med.df, aes(xintercept = dist_norm_rbls), lty = 4) +
  facet_wrap(~species_to, scales = 'free') +
  xlim(0, 5) +
  labs(title = 'A to sp. seq to seq normalised distance')
# dev.off()

A_sort <- med.df[order(med.df$dist_norm_mrca), 'species_to']

sp_dat <- read.csv('../outputs/sptree_dist.csv')
a <- med.df[, c(2, 3, 5, 7, 9, 10, 12)] - sp_dat[1:4, ]$dist
apply(a[-4, ], 2, sum)

write.csv(a, file = '../outputs/norm_tau_differences.csv')
