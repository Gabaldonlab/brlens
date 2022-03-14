library(ggplot2)
library(ggpubr)

theme_set(theme_bw())

phylomes <- c('ECOLI' = '0677', 'YEAST' = '0907', 'HUMAN' = '0739',
              'YEAST' = '0005', 'HUMAN' = '0076')

plots <- c()
for (i in 1:length(phylomes)) {
  phylome <- phylomes[i]
  refsp <- names(phylomes)[i]
  
  dat <- read.csv(paste0('../../11_cladenorm/outputs/', phylome, '_norm.csv'))
  
  p <- ggplot(dat, aes(mrca_to_tip_median, col = Normalising.group)) +
    geom_density() +
    xlim(0, quantile(dat$mrca_to_tip_median, 0.99)) +
    xlab('Group MRCA to tip median') +
    labs(title = paste(refsp, 'phylome'))
  
  plots[[phylome]] <- p
  
  pdf(paste0('../outputs/', phylome, '_normfact_dens.pdf'),
      width = 6.5, height = 3.5)
  print(p)
  dev.off()
}

pdf('../outputs/qfo_normfact_dens.pdf', width = 10, height = 3)
ggarrange(plots[['0677']], plots[['0907']], plots[['0739']], labels = 'auto',
          align = 'hv', common.legend = TRUE, ncol = 3)
dev.off()
