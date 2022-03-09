# Analysis of phylome distances normalised by clade
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
phylome <- '0005'
refsp <- 'YEAST'

# Loading phylome info
phyinfo <- read.csv(paste0('../../11_cladenorm/data/', phylome, '_norm_groups.csv'))

# Loading and basic data management of phylome distances
dat <- read.csv(paste0('../../11_cladenorm/outputs/', phylome, '_dist.csv'))
spdat <- dat[which(dat$mrca_type == 'S' & dat$from_sp == refsp | 
                     (dat$to_sp == refsp & dat$from_sp != dat$to_sp)), ]

spdat$sp_to <- apply(spdat, 1, get_other, ref = refsp)

sptreedat <- read.csv(paste0('../../11_cladenorm/outputs/', phylome,
                             '_sptree_dist.csv'))
sptreedat$sp_to <- apply(sptreedat, 1, get_other, ref = refsp)

# Plotting data ----
# Descriptive plots - densities
a <- ggplot(spdat, aes(x = dist, colour = sp_to)) +
  geom_density(show.legend = FALSE) +
  xlim(0, 5) +
  xlab('Raw distance')

b <- ggplot(spdat, aes(x = ndist_A, colour = sp_to)) +
  geom_density(show.legend = FALSE) +
  xlim(0, 7.5) +
  xlab('Clade A normalised distance')

ggarrange(a, b, labels = 'auto')

ggplot(spdat, aes(x = dist, colour = sp_to, fill = sp_to)) +
  geom_density(show.legend = FALSE, alpha = 0.6) +
  facet_wrap(~sp_to, scales = 'free') +
  labs(title = paste(refsp, 'to sp raw distance'))

ggplot(spdat, aes(x = ndist_A, colour = sp_to, fill = sp_to)) +
  geom_density(show.legend = FALSE, alpha = 0.6) +
  facet_wrap(~sp_to, scales = 'free') +
  labs(title = paste(refsp, 'to sp - clade A normalised distance'))

phymeds <- spdat %>%
  group_by(sp_to) %>%
  summarise('dist' = median(dist),
            'ndist_A' = median(ndist_A))

sptd <- sptreedat[which((sptreedat$from_sp == refsp | sptreedat$to_sp == refsp) &
                          sptreedat$from_sp != sptreedat$to_sp), ]
a <- merge(phymeds, sptd, by = 'sp_to', suffixes = c('_phy', '_sp'))

ggplot(a, aes(dist_sp, ndist_A_phy)) +
  geom_point()
