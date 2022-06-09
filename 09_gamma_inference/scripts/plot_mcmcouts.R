# Plot MCMC of all the posteriors
# Moisès Bernabeu
# Barcelona, May 2022

library(ggplot2)
library(ggpubr)

theme_set(theme_bw())

# Import parsed data
load('../outputs/sp2sp_mcmc_gathered.RData')

# humcmcg <- humcmcg[-which(humcmcg$param == 'mo' & humcmcg$value == 0), ]
humcmcg <- humcmcg[-which(humcmcg$spto %in% c('PANTR', 'MOUSE', 'MONDO')), ]

# Plotting posterior densities
# pdf('../outputs/human_allpars_dens.pdf', width = 8, height = 4)
ggplot(humcmcg, aes(x = value, color = spto)) +
  geom_density() +
  facet_wrap(~param, scales = 'free') +
  # ylim(0, 100) +
  theme(legend.position = 'none')
  # guides(color = guide_legend(nrow = 4, byrow = TRUE))
# dev.off()

# pdf('../outputs/yeast_allpars_dens.pdf', width = 8, height = 4)
ggplot(yemcmcg, aes(x = value, color = spto)) +
  geom_density() +
  facet_wrap(~param, scales = 'free') +
  # ylim(0, 100) +
  theme(legend.position = 'none')
  # guides(color = guide_legend(nrow = 2, byrow = TRUE))
# dev.off()

# Mode posterior densities
ggplot(humcmcg[which(humcmcg$param == 'mo'), ], aes(x = value, color = spto)) +
  geom_density()

a <- ggplot(humcmcg[which(humcmcg$param == 'mo'), ], aes(x = value, color = spto)) +
  geom_density() +
  xlab('Mode') 

b <- ggplot(humcmcg[which(humcmcg$param == 'mo'), ], aes(x = value, color = spto)) +
  geom_density() +
  xlim(0, 0.2) +
  # ylim(0, 100) +
  xlab('Mode') +
  ylab('Density')

c <- ggplot(humcmcg[which(humcmcg$param == 'mo'), ], aes(x = value, color = spto)) +
  geom_density() +
  xlim(1.5, 5) +
  xlab('Mode') +
  ylab('Density')

# pdf('../outputs/human_mode_postdens.pdf', width = 10.8, height = 4.5)
ggarrange(a, b, c, common.legend = TRUE, legend = 'none',
          ncol = 3, align = 'hv', labels = 'auto')
# dev.off()

# pdf('../outputs/human_mode_postdens.pdf', width = 7.5, height = 2.7)
ggarrange(b, c, common.legend = TRUE, legend = 'none',
          ncol = 2, align = 'hv', labels = 'auto')
# dev.off()

a <- ggplot(yemcmcg[which(yemcmcg$param == 'mo'), ], aes(x = value, color = spto)) +
  geom_density() +
  xlab('Mode')

b <- ggplot(yemcmcg[which(yemcmcg$param == 'mo'), ], aes(x = value, color = spto)) +
  geom_density() +
  xlim(0, 0.25) +
  xlab('Mode') +
  ylab('Density')

c <- ggplot(yemcmcg[which(yemcmcg$param == 'mo'), ], aes(x = value, color = spto)) +
  geom_density() +
  xlim(1, 3.5) +
  xlab('Mode') +
  ylab('Density')

# pdf('../outputs/yeast_mode_postdens.pdf', width = 10.8, height = 3.65)
ggarrange(a, b, c, common.legend = TRUE, legend = 'bottom',
          ncol = 3, align = 'hv', labels = 'auto')
# dev.off()

# pdf('../outputs/yeast_mode_postdens.pdf', width = 7.5, height = 2.7)
ggarrange(b, c, common.legend = TRUE, legend = 'none',
          ncol = 2, align = 'hv', labels = 'auto')
# dev.off()

load('../outputs/event_mcmc_gathered.RData')

# pdf('../outputs/event_allpars_dens.pdf', width = 7.5, height = 4.40)
ggplot(evmcmcg, aes(x = value, color = spto)) +
  geom_density() +
  facet_wrap(~param, scales = 'free') +
  theme(legend.position = 'bottom') +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
# dev.off()
