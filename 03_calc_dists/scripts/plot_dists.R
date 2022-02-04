# plot_dists.R -- Useful plots for distances
# 
# The program gets the calc_dists.py csv output and plot densities, histograms
# and boxplots for the relation between duplications, speciations and
# distances.
# 
# Requirements: ggplot2
# 
# Written by Moisès Bernabeu moigil.bernabeu@gmail.com
# January 2022

# Loading packages ----
library(e1071)
library(ggplot2)
library(GGally)
# library(ggpubr)
library(patchwork)

# Settings ----
theme_set(theme_bw())

# Loading data ----
files <- list.files('../outputs', pattern = 'all', full.names = TRUE)

# Functions definition ----
sumstats <- function(x) {
  c(summary(x), kurt = kurtosis(x, na.rm = TRUE),
    skew = skewness(x, na.rm = TRUE))
}

sumdf <- rbind(sumstats(rnorm(200)))
sumdf <- sumdf[-1, ]

for (file in files) {
  print(file)
  phylome <- strsplit(tail(strsplit(file, '/')[[1]], n = 1), '_')[[1]][1]
  dat <- read.csv(file)
  
  # Basic stats ----
  sum <- rbind(sumstats(dat$seed_dist),
               sumstats(dat$seed_sp),
               sumstats(dat$seed_dupl),
               sumstats(dat$og_dist),
               sumstats(dat$og_sp),
               sumstats(dat$og_dupl))
  row.names(sum) <- c(paste(phylome, 'seed_dist', sep = '_'),
                      paste(phylome, 'seed_sp', sep = '_'),
                      paste(phylome, 'seed_dupl', sep = '_'),
                      paste(phylome, 'og_dist', sep = '_'),
                      paste(phylome, 'og_sp', sep = '_'),
                      paste(phylome, 'og_dupl', sep = '_'))
  sumdf <- rbind(sumdf, sum)
  
  # Histograms ----
  # Seed to leaf speciation events
  slsp <- ggplot(dat, aes(seed_sp)) +
    geom_bar(col = 'black', fill = 'black', alpha = 0.4) +
    xlab('Seed-leaf speciations')
  
  # Seed to leaf duplication events
  sldp <- ggplot(dat, aes(seed_dupl)) +
    geom_bar(col = 'black', fill = 'black', alpha = 0.4) +
    xlab('Seed-leaf duplications')
  
  # Outgroup to leaf speciation events
  olsp <- ggplot(dat, aes(og_sp)) +
    geom_bar(col = 'black', fill = 'black', alpha = 0.4) +
    xlab('Outgroup-leaf speciations')
  
  # Outgroup to leaf duplication events
  oldp <- ggplot(dat, aes(og_dupl)) +
    geom_bar(col = 'black', fill = 'black', alpha = 0.4) +
    xlab('Outgroup-leaf duplications')
  
  # Joining the plots
  pdf(paste0('../outputs/', phylome, '_sp_dp_hist.pdf'),
      width = 7, height = 4.5)
  print((slsp + sldp) / (olsp + oldp))
  # print(slsp)
  # print(sldp)
  # print(olsp)
  # print(oldp)
  dev.off()
  
  # Density of numeric variables ----
  sd <- ggplot(dat, aes(seed_dist)) +
    geom_density(col = 'black', fill = 'black', alpha = 0.4) +
    xlab('Seed to leaf distance') +
    xlim(0, 30)
  
  od <- ggplot(dat, aes(og_dist)) +
    geom_density(col = 'black', fill = 'black', alpha = 0.4) +
    xlab('Outgroup to leaf distance') +
    xlim(0, 30)
  
  pdf(paste0('../outputs/', phylome, '_dens.pdf'), width = 7, height = 2.25)
  print(sd + od)
  # print(sd)
  # print(od)
  dev.off()
  
  # pdf(paste0('../outputs/', phylome, '_pairs.pdf'), width = 12, height = 12)
  # print(ggpairs(dat[, 5:10]))
  # dev.off()
}

write.csv(sumdf, file = '../outputs/sumstats.csv')
