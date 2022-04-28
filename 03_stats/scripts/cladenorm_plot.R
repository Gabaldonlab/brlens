# Analysis of phylome distances normalised by clade
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
phylomes <- c('YEAST' = '0005', 'HUMAN' = '0076')

for (i in 1:length(phylomes)) {
  phylome <- phylomes[i]
  refsp <- names(phylomes)[i]
  
  # Loading phylome info
  phyinfo <- read.csv(paste0('../../02_get_distances/data/', phylome, '_norm_groups.csv'))
  
  # Loading and basic data management of phylome distances
  dat <- read.csv(paste0('../../02_get_distances/outputs/', phylome, '_dist.csv'))
  spdat <- dat[which(dat$mrca_type == 'S' & dat$from_sp == refsp | 
                       (dat$to_sp == refsp & dat$from_sp != dat$to_sp)), ]
  
  spdat$sp_to <- apply(spdat, 1, get_other, ref = refsp)
  
  sptreedat <- read.csv(paste0('../../02_get_distances/outputs/', phylome,
                               '_sptree_dist.csv'))
  sptreedat$sp_to <- apply(sptreedat, 1, get_other, ref = refsp)
  
  # Plotting data ----
  # Descriptive plots - densities
  a <- ggplot(spdat, aes(x = dist, colour = sp_to)) +
    geom_density(show.legend = FALSE) +
    xlim(0, quantile(spdat$dist, 0.95)) +
    xlab('Raw distance')
  
  b <- ggplot(spdat, aes(x = ndist_A, colour = sp_to)) +
    geom_density(show.legend = FALSE) +
    xlim(0, quantile(spdat$ndist_A, 0.9)) +
    xlab('Clade A normalised distance')
  
  pdf(paste0('../outputs/', phylome, '_densities.pdf'), width = 9, height = 3)
  print(ggarrange(a, b, labels = 'auto'))
  dev.off()
  
  # ggplot(spdat, aes(x = dist, colour = sp_to, fill = sp_to)) +
  #   geom_density(show.legend = FALSE, alpha = 0.6) +
  #   xlim(0, quantile(spdat$dist, 0.95)) +
  #   facet_wrap(~sp_to, scales = 'free') +
  #   labs(title = paste(refsp, 'to sp raw distance'))
  # 
  # ggplot(spdat, aes(x = ndist_A, colour = sp_to, fill = sp_to)) +
  #   geom_density(show.legend = FALSE, alpha = 0.6) +
  #   xlim(0, quantile(spdat$ndist_A, 0.9)) +
  #   facet_wrap(~sp_to, scales = 'free') +
  #   labs(title = paste(refsp, 'to sp - clade A normalised distance'))
  
  phymeds <- spdat %>%
    group_by(sp_to) %>%
    summarise('dist' = median(dist),
              'ndist_A' = median(ndist_A))
  
  sptd <- sptreedat[which((sptreedat$from_sp == refsp | sptreedat$to_sp == refsp) &
                            sptreedat$from_sp != sptreedat$to_sp), ]
  spphydf <- merge(phymeds, sptd, by = 'sp_to', suffixes = c('_phy', '_sp'))
  
  compp <- ggplot(spphydf, aes(dist_sp, ndist_A_phy)) +
    geom_point() +
    xlab('Species tree distance') +
    ylab('Phylome norm. distance median')
  compp
  
  # All comparisons ----
  comp_dat <- dat[which(dat$dupl == 0 & dat$mrca_type == 'S'), ]
  comp_dat$from_to <- paste(comp_dat$from_sp, comp_dat$to_sp, sep = '-')
  comp_sptree <- sptreedat
  comp_sptree$from_to <- paste(comp_sptree$from_sp, comp_sptree$to_sp, sep = '-')
  
  comp_phy <- comp_dat %>%
    group_by(from_to) %>%
    summarise('dist' = median(dist),
              'ndist_A' = median(ndist_A, na.rm = TRUE),
              'ndist_B' = median(ndist_B, na.rm = TRUE),
              'sp' = mean(sp),
              'dupl' = mean(dupl))
  
  comp <- merge(comp_phy, comp_sptree, by = 'from_to',
                suffixes = c('_phy', '_sp'))
  
  compp_all <- ggplot(comp, aes(dist_sp, ndist_A_phy)) +
    geom_point() +
    xlab('Species tree distance') +
    ylab('Phylome norm. distance median')
  compp_all
  
  pdf(paste0('../outputs/', phylome, '_phy_vs_sptree.pdf'), width = 9, height = 3)
  print(ggarrange(compp, compp_all, align = 'hv', labels = 'auto'))
  dev.off()
  
  # Duplications and speciations ----
  dpsp_dat <- dat[which(dat$mrca_type == 'S'), ]
  dpsp_dat$from_to <- paste(dpsp_dat$from_sp, dpsp_dat$to_sp, sep = '-')
  
  dpsp_dat <- dpsp_dat %>%
    group_by(from_to) %>%
    summarise('dist' = median(dist),
              'ndist_A' = median(ndist_A, na.rm = TRUE),
              'ndist_B' = median(ndist_B, na.rm = TRUE),
              'sp' = mean(sp),
              'dupl' = mean(dupl))
  
  # Raw distance
  a <- ggplot(dpsp_dat, aes(dist, dupl)) +
    geom_point() +
    xlab('Pairwise distance') +
    ylab('Number of duplication events')
  
  b <- ggplot(dpsp_dat, aes(dist, sp)) +
    geom_point() +
    xlab('Pairwise distance') +
    ylab('Number of speciation events')
  
  pdf(paste0('../outputs/', phylome, '_dist_dupl_sp_.pdf'), width = 9, height = 3)
  print(ggarrange(a, b, align = 'hv', labels = 'auto'))
  dev.off()
  
  # Normalised distance
  a <- ggplot(dpsp_dat, aes(ndist_A, dupl)) +
    geom_point() +
    xlab('Clade A normalised pairwise distance') +
    ylab('Number of duplication events')
  
  b <- ggplot(dpsp_dat, aes(ndist_A, sp)) +
    geom_point() +
    xlab('Clade A normalised pairwise distance') +
    ylab('Number of speciation events')
  
  pdf(paste0('../outputs/', phylome, '_ndist_dupl_sp_.pdf'), width = 9, height = 3)
  print(ggarrange(a, b, align = 'hv', labels = 'auto'))
  dev.off()
}

# Special plots ----
yedat <- read.csv('../../02_get_distances/outputs/0005_dist.csv')
yedat <- yedat[which(yedat$mrca_type == 'S' & (yedat$from_sp == 'YEAST' | 
                     yedat$to_sp == 'YEAST') & yedat$from_sp != yedat$to_sp), ]
yedat$sp_to <- apply(yedat, 1, get_other, ref = 'YEAST')

str(yedat)
summary(yedat)

hudat <- read.csv('../../02_get_distances/outputs/0076_dist.csv')
hudat <- hudat[which(hudat$mrca_type == 'S' & (hudat$from_sp == 'HUMAN' | 
                       hudat$to_sp == 'HUMAN') & hudat$from_sp != hudat$to_sp), ]
hudat$sp_to <- apply(hudat, 1, get_other, ref = 'HUMAN')

str(hudat)
summary(hudat)

yerdens <- ggplot(yedat, aes(dist, col = sp_to)) +
  geom_density(show.legend = FALSE) +
  xlim(0, quantile(yedat$dist, 0.99)) +
  xlab('Yeast to sp. distance') +
  ylab('Density')

yerdens

hurdens <- ggplot(hudat, aes(dist, col = sp_to)) +
  geom_density(show.legend = FALSE) +
  xlim(0, quantile(hudat$dist, 0.99)) +
  xlab('Human to sp. distance') +
  ylab('Density')

hurdens

# pdf('../outputs/raw_densities.pdf', width = 6.3, height = 2.3)
ggarrange(yerdens, hurdens, align = 'hv', labels = 'auto')
# dev.off()

yendens <- ggplot(yedat, aes(ndist_A, col = sp_to)) +
  geom_density(show.legend = FALSE) +
  xlim(0, quantile(yedat$dist, 0.999)) +
  xlab('Yeast to sp. norm. distance') +
  ylab('Density')

yendens

hundens <- ggplot(hudat, aes(ndist_A, col = sp_to)) +
  geom_density(show.legend = FALSE) +
  xlim(0, quantile(hudat$dist, 0.99)) +
  xlab('Human to sp. norm. distance') +
  ylab('Density')

hundens

# pdf('../outputs/norm_densities.pdf', width = 6.3, height = 2.3)
ggarrange(yendens, hundens, align = 'hv', labels = 'auto')
# dev.off()

