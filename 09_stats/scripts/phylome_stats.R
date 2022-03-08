# Analysis of phylome distances
# Moisès Bernabeu
# March, 2022

# Loading libraries ----
library(ggplot2)
library(dplyr)

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

median_sp <- function(x, df) {
  by(x, df[, 'sp_to'], median, na.rm = TRUE)
}

# Import data ----
phylome <- '0739'
refsp <- 'HUMAN'

# Loading phylome info
phyinfo <- read.csv('../../07_normalisation/data/README',
                    sep = '\t', row.names = 3)

# Loading and basic data management of phylome distances
dat <- read.csv(paste0('../data/', phylome, '_dist.csv'))
spdat <- dat[which(dat$mrca_type == 'S' & dat$from_sp == refsp | 
                     dat$to_sp == refsp & dat$from_sp != dat$to_sp), ]

spdat$sp_to <- apply(spdat, 1, get_other, ref = refsp)
spdat$kingdom_to <- phyinfo[spdat$sp_to, 'kingdom']

# Loading normalisation parameters
nor <- read.csv(paste0('../data/', phylome, '_norm.csv'))

# Plotting data ----
# Descriptive plots - densities
ggplot(spdat, aes(x = dist, colour = sp_to)) +
  geom_density(show.legend = FALSE) +
  xlim(0, 15)

ggplot(spdat, aes(x = mrca_ndist, colour = sp_to)) +
  geom_density(show.legend = FALSE) +
  xlim(0, 10)

ggplot(spdat, aes(x = rbl_ndist, colour = sp_to)) +
  geom_density(show.legend = FALSE) +
  xlim(0, 10)

ggplot(spdat, aes(x = dist, colour = sp_to, fill = sp_to)) +
  geom_density(show.legend = FALSE, alpha = 0.6) +
  facet_wrap(~sp_to, scales = 'free') +
  labs(title = paste(refsp, 'to sp - QfO phylome (78 species)'))

ggplot(spdat, aes(x = mrca_ndist, colour = sp_to, fill = sp_to)) +
  geom_density(show.legend = FALSE, alpha = 0.6) +
  facet_wrap(~sp_to, scales = 'free') +
  labs(title = paste(refsp, 'to sp - QfO phylome (78 species)'))

ggplot(spdat, aes(x = rbl_ndist, colour = sp_to, fill = sp_to)) +
  geom_density(show.legend = FALSE, alpha = 0.6) +
  facet_wrap(~sp_to, scales = 'free') +
  labs(title = paste(refsp, 'to sp - QfO phylome (78 species)'))

spdat[, c('tree', 'dist', 'mrca_ndist', 'rbl_ndist')] %>%
  group_by(tree) %>%
  summarise(median(dist), median(mrca_dist))
