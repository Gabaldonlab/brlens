# Analysis of phylome distances normalised by clade - QfO
# Moisès Bernabeu
# March, 2022

# Loading libraries ----
library(ggplot2)
library(dplyr)
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

median_sp <- function(x, df) {
  by(x, df[, 'sp_to'], median, na.rm = TRUE)
}

# Import data ----
phylomes <- c('ECOLI' = '0677', 'YEAST' = 'O907', 'HUMAN' = '0739')

phylome <- '0677'
refsp <- 'ECOLI'

# Loading phylome info
phyinfo <- read.csv(paste0('../../11_cladenorm/data/qfo_78_norm_groups.csv'))

# Loading and basic data management of phylome distances
dat <- read.csv(paste0('../../11_cladenorm/outputs/', phylome, '_dist.csv'))
spdat <- dat[which(dat$mrca_type == 'S' & dat$from_sp == refsp | 
                     (dat$to_sp == refsp & dat$from_sp != dat$to_sp)), ]

spdat$sp_to <- apply(spdat, 1, get_other, ref = refsp)

sptreedat <- read.csv(paste0('../../11_cladenorm/outputs/qfo78_sp_tree_dist.csv'))
sptreedat$sp_to <- apply(sptreedat, 1, get_other, ref = refsp)

# Plotting data ----
# Descriptive plots - densities
a <- ggplot(spdat, aes(x = dist, colour = sp_to)) +
  geom_density(show.legend = FALSE) +
  xlim(0, 25) +
  xlab('Raw distance')

b <- ggplot(spdat, aes(x = ndist_eukaryota, colour = sp_to)) +
  geom_density(show.legend = FALSE) +
  xlim(0, 30) +
  xlab('Eukaryota clade normalised distance')

ggarrange(a, b, labels = 'auto')

a <- ggplot(spdat, aes(x = ndist_eukaryota, colour = sp_to)) +
  geom_density(show.legend = FALSE) +
  xlim(0, 30) +
  xlab('Eukaryota clade normalised distance') +
  labs(title = paste(refsp, 'phylome'))
a

b <- ggplot(spdat, aes(x = ndist_bacteria, colour = sp_to)) +
  geom_density(show.legend = FALSE) +
  xlim(0, 30) +
  xlab('Bacteria clade normalised distance') +
  labs(title = paste(refsp, 'phylome'))
b

c <- ggplot(spdat, aes(x = ndist_archaea, colour = sp_to)) +
  geom_density(show.legend = FALSE) +
  xlim(0, 30) +
  xlab('Archaea clade normalised distance') +
  labs(title = paste(refsp, 'phylome'))
c

ggarrange(a, b, c, align = 'h', nrow = 1)

phymeds <- spdat %>%
  group_by(sp_to) %>%
  summarise('dist' = median(dist),
            'ndist_eukaryota' = median(ndist_eukaryota, na.rm = TRUE))

sptd <- sptreedat[which((sptreedat$from_sp == refsp | sptreedat$to_sp == refsp) &
                          sptreedat$from_sp != sptreedat$to_sp), ]
a <- merge(phymeds, sptd, by = 'sp_to', suffixes = c('_phy', '_sp'))
row.names(a) <- a$sp_to
row.names(phyinfo) <- phyinfo$Proteome

a$domain <- phyinfo[row.names(a), 'kingdom']

compp <- ggplot(a, aes(dist_sp, ndist_eukaryota_phy, colour = domain)) +
  geom_point()
compp

# All comparisons ----
comp_dat <- dat[which(dat$dupl == 0 & dat$mrca_type == 'S'), ]
comp_dat$from_to <- paste(comp_dat$from_sp, comp_dat$to_sp, sep = '-')
comp_sptree <- sptreedat
comp_sptree$from_to <- paste(comp_sptree$from_sp, comp_sptree$to_sp, sep = '-')

comp_phy <- comp_dat %>%
  group_by(from_to) %>%
  summarise('dist' = median(dist),
            'ndist_eukaryota' = median(ndist_eukaryota, na.rm = TRUE),
            'ndist_bacteria' = median(ndist_bacteria, na.rm = TRUE),
            'ndist_archaea' = median(ndist_archaea, na.rm = TRUE),
            'sp' = mean(sp),
            'dupl' = mean(dupl))

comp <- merge(comp_phy, comp_sptree, by = 'from_to',
              suffixes = c('_phy', '_sp'))

compp_all <- ggplot(comp, aes(dist_sp, ndist_eukaryota_phy)) +
  geom_point()

ggarrange(compp, compp_all, align = 'hv', labels = 'auto', common.legend = TRUE)
