# Data analysis of phylome distances - Normalisation
# Moisès Bernabeu
# Barcelona, February, 2022

# Libraries loading ----
library(ggplot2)
library(ggpubr)
library(tidyr)

theme_set(theme_bw())

# Functions definition ----
get_other <- function(x, ref) {
  if (x['from_sp'] != ref) {
    y <- x['from_sp']
  } else {
    y <- x['to_sp']
  }
  return(y)
}

median_sp <- function(x, df) {
  by(x, df[, 'sp_to'], median, na.rm = TRUE)
}

# Normalised distances dataframe generation ----
phylome <- '0076'
refsp <- 'HUMAN'

dat <- read.csv(paste0('../data/', phylome, '_dist_noh.csv'))
spdat <- dat[which(dat$mrca_type == 'S' & dat$from_sp == refsp | 
                     dat$to_sp == refsp & dat$from_sp != dat$to_sp), ]
meds <- apply(spdat[!duplicated(spdat$tree), c(11:16)], 2, median)

oldvarnms <- names(spdat)[-15][11:15]
newvarnms <- c('st_ndist', 'mrca_ndist', 'root_ndist',
               'twdth_ndist', 'brl_ndist')

distdf <- data.frame('sp_to' = apply(spdat, 1, get_other, ref = refsp))

distdf$raw_dist <- spdat$dist
for (i in 1:length(oldvarnms)) {
  distdf[, newvarnms[i]] <- spdat$dist / (spdat[, oldvarnms[i]] /
                                            meds[oldvarnms[i]])
}

distdfg <- gather(distdf, key = 'distance', value = 'value', -sp_to)

# Importing species tree data ----
spref <- read.csv(paste0('../data/', phylome, '_sptree_sptree_dist.csv'))
spref <- spref[which(spref$from_sp == refsp | spref$to_sp == refsp), ]
spref$sp_to <- apply(spref, 1, get_other, ref = refsp)

spdists <- data.frame('sp_to' = spref$sp_to, 'dist_sp' = spref$dist)
phydists <- data.frame(apply(distdf[, -1], 2,
                             FUN = median_sp, df = distdf)[spdists$sp_to, ])

sist_group <- read.csv(paste0('../data/', phylome, '_sister_group.csv'),
                       row.names = 1)

sp_vs_phy <- cbind(spdists, phydists,
                   'group' = factor(sist_group[row.names(phydists), ]))

# Writing files ----
write.csv(distdf, file = paste0('../outputs/', phylome, '_distdf.csv'))
write.csv(sp_vs_phy, file = paste0('../outputs/', phylome, '_sp_vs_phylome.csv'))

# Plots ----
distdf <- read.csv(paste0('../outputs/', phylome, '_distdf.csv'),
                   row.names = 1)
sp_vs_phy <- read.csv(paste0('../outputs/', phylome, '_sp_vs_phylome.csv'),
                      row.names = 1)
sp_vs_phy$group <- factor(sp_vs_phy$group)

pdf(paste0('../outputs/', phylome, '_all_dens.pdf'),
    width = 6, height = 4)
ggplot(distdf, aes(x = raw_dist, colour = sp_to)) +
  geom_density(show.legend = FALSE) +
  xlim(0, 5) +
  labs(title = paste0(refsp, ' to species'))
dev.off()

pdf(paste0('../outputs/', phylome, '_all_norm_dens.pdf'),
    width = 6, height = 4)
ggplot(distdf, aes(x = brl_ndist, colour = sp_to)) +
  geom_density(show.legend = FALSE) +
  xlim(0, 7) +
  labs(title = paste0(refsp, ' to species'))
dev.off()

pdf(paste0('../outputs/', phylome, '_sp_vs_phylome.pdf'),
    width = 5.5, height = 5)
for (varnm in names(sp_vs_phy)[-c(1:2, 9)]) {
  p <- ggplot(sp_vs_phy, aes(x = dist_sp, y = get(varnm), colour = group)) +
    geom_point() +
    ylab(varnm) +
    theme(legend.position = 'bottom') +
    labs(title = paste0(refsp, ' to species'))
  print(p)
}
dev.off()

pdf(paste0('../outputs/', phylome, '_dist_dens.pdf'),
    width = 10, height = 6)
for (varnm in names(distdf)[-1]) {
  p <- ggplot(distdf, aes(x = get(varnm), colour = sp_to, fill = sp_to)) +
    geom_density(show.legend = FALSE, alpha = 0.6) +
    xlab(varnm) +
    xlim(0, quantile(x = distdf[, varnm], probs = 0.95, na.rm = TRUE)) +
    facet_wrap(~sp_to, scales = 'free')
  print(p)
}
dev.off()
