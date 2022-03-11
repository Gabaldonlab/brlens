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
phylomes <- c('ECOLI' = '0677', 'YEAST' = '0907', 'HUMAN' = '0739')

# Loading phylome info
phyinfo <- read.csv(paste0('../../11_cladenorm/data/qfo_78_norm_groups.csv'))

sptreedat <- read.csv(paste0('../../11_cladenorm/outputs/qfo78_sp_tree_dist.csv'))

for (i in 1:length(phylomes)) {
  phylome <- phylomes[i]
  refsp <- names(phylomes)[i]
  
  # Loading and basic data management of phylome distances
  dat <- read.csv(paste0('../../11_cladenorm/outputs/', phylome, '_dist.csv'))
  spdat <- dat[which(dat$mrca_type == 'S' & dat$from_sp == refsp | 
                       (dat$to_sp == refsp & dat$from_sp != dat$to_sp)), ]
  
  spdat$sp_to <- apply(spdat, 1, get_other, ref = refsp)
  sptreedat$sp_to <- apply(sptreedat, 1, get_other, ref = refsp)
  
  # Plotting data ----
  # Descriptive plots - densities
  a <- ggplot(spdat, aes(x = dist, colour = sp_to)) +
    geom_density(show.legend = FALSE) +
    xlim(0, quantile(spdat$dist, 0.95)) +
    xlab('Eukaryota clade normalised distance') +
    labs(title = paste(refsp, 'phylome'))
  
  b <- ggplot(spdat, aes(x = ndist_bacteria, colour = sp_to)) +
    geom_density(show.legend = FALSE) +
    xlim(0, quantile(spdat$ndist_bacteria, 0.95, na.rm = TRUE)) +
    xlab('Bacteria clade normalised distance') +
    labs(title = paste(refsp, 'phylome'))

  c <- ggplot(spdat, aes(x = ndist_archaea, colour = sp_to)) +
    geom_density(show.legend = FALSE) +
    xlim(0, quantile(spdat$ndist_archaea, 0.95, na.rm = TRUE)) +
    xlab('Archaea clade normalised distance') +
    labs(title = paste(refsp, 'phylome'))
  
  d <- ggplot(spdat, aes(x = ndist_eukaryota, colour = sp_to)) +
    geom_density(show.legend = FALSE) +
    xlim(0, quantile(spdat$ndist_eukaryota, 0.9, na.rm = TRUE)) +
    xlab('Eukaryota clade normalised distance') +
    labs(title = paste(refsp, 'phylome'))

  pdf(paste0('../outputs/', phylome, '_densities.pdf'), width = 8, height = 6)
  print(ggarrange(a, b, c, d, align = 'hv'))
  dev.off()
  
  phymeds <- spdat %>%
    group_by(sp_to) %>%
    summarise('dist' = median(dist),
              'ndist_eukaryota' = median(ndist_eukaryota, na.rm = TRUE),
              'ndist_bacteria' = median(ndist_bacteria, na.rm = TRUE),
              'ndist_archaea' = median(ndist_archaea, na.rm = TRUE),
              'sp' = mean(sp),
              'dupl' = mean(dupl))
  
  sptd <- sptreedat[which((sptreedat$from_sp == refsp | sptreedat$to_sp == refsp) &
                            sptreedat$from_sp != sptreedat$to_sp), ]
  spphydf <- merge(phymeds, sptd, by = 'sp_to', suffixes = c('_phy', '_sp'))
  row.names(spphydf) <- spphydf$sp_to
  row.names(phyinfo) <- phyinfo$Proteome
  
  spphydf$domain <- phyinfo[row.names(spphydf), 'kingdom']
  
  a <- ggplot(spphydf, aes(dist_sp, ndist_eukaryota_phy, colour = domain)) +
    geom_point() +
    xlab('Species tree distance') +
    ylab('Phylome euk. norm. distance med.')

  b <- ggplot(spphydf, aes(dist_sp, ndist_bacteria_phy, colour = domain)) +
    geom_point() +
    xlab('Species tree distance') +
    ylab('Phylome bac. norm. distance med.')

  c <- ggplot(spphydf, aes(dist_sp, ndist_archaea_phy, colour = domain)) +
    geom_point() +
    xlab('Species tree distance') +
    ylab('Phylome arch. norm. distance med.')
  
  pdf(paste0('../outputs/', phylome, '_seed_to_sp_phy_sptree.pdf'), width = 15, height = 3.5)
  print(ggarrange(a, b, c, align = 'hv', labels = 'auto', nrow = 1, common.legend = TRUE))
  dev.off()
  
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
  
  a <- ggplot(comp, aes(dist_sp, ndist_eukaryota_phy)) +
    geom_point() +
    xlab('Species tree distance') +
    ylab('Phylome euk. norm. distance med.')
  
  b <- ggplot(comp, aes(dist_sp, ndist_bacteria_phy)) +
    geom_point() +
    xlab('Species tree distance') +
    ylab('Phylome bac. norm. distance med.')
  
  c <- ggplot(comp, aes(dist_sp, ndist_archaea_phy)) +
    geom_point() +
    xlab('Species tree distance') +
    ylab('Phylome arch. norm. distance med.')
  
  pdf(paste0('../outputs/', phylome, '_pw_phy_sptree.pdf'), width = 15, height = 3.5)
  print(ggarrange(a, b, c, align = 'hv', labels = 'auto', nrow = 1, common.legend = TRUE))
  dev.off()
  
  comp_dat <- dat[which(dat$mrca_type == 'S'), ]
  comp_dat$from_to <- paste(comp_dat$from_sp, comp_dat$to_sp, sep = '-')
  comp_phy <- comp_dat %>%
    group_by(from_to) %>%
    summarise('dist' = median(dist),
              'ndist_eukaryota' = median(ndist_eukaryota, na.rm = TRUE),
              'ndist_bacteria' = median(ndist_bacteria, na.rm = TRUE),
              'ndist_archaea' = median(ndist_archaea, na.rm = TRUE),
              'sp' = mean(sp),
              'dupl' = mean(dupl))
  
  a <- ggplot(comp_phy, aes(ndist_bacteria, dupl)) +
    geom_point() +
    xlab('Bac. normalised pairwise distance') +
    ylab('Number of duplication events')
  
  b <- ggplot(comp_phy, aes(ndist_archaea, dupl)) +
    geom_point() +
    xlab('Arc. normalised pairwise distance') +
    ylab('Number of duplication events')
  
  c <- ggplot(comp_phy, aes(ndist_eukaryota, dupl)) +
    geom_point() +
    xlab('Euk. normalised pairwise distance') +
    ylab('Number of duplication events')

  d <- ggplot(comp_phy, aes(ndist_bacteria, sp)) +
    geom_point() +
    xlab('Bac. normalised pairwise distance') +
    ylab('Number of speciation events')
  
  e <- ggplot(comp_phy, aes(ndist_archaea, sp)) +
    geom_point() +
    xlab('Arc. normalised pairwise distance') +
    ylab('Number of speciation events')
  
  f <- ggplot(comp_phy, aes(ndist_eukaryota, sp)) +
    geom_point() +
    xlab('Euk. normalised pairwise distance') +
    ylab('Number of speciation events')
    
  pdf(paste0('../outputs/', phylome, '_sp_dupl.pdf'), width = 15, height = 7)
  print(ggarrange(a, b, c, d, e, f, align = 'hv', labels = 'auto',
                  nrow = 2, ncol = 3, common.legend = TRUE))
  dev.off()
}
