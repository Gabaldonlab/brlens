# Plottin inferred Gamma distributions
# Moisès Bernabeu
# Barcelona, June 2022

# Loading packages ----
library(ggplot2)

theme_set(theme_bw())

# Defining functions ----
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Loading data ----
load('../data/event_dist.RData')
load('../outputs/event_mcmc_gathered.RData')

evpars <- by(evmcmcg$value, paste(evmcmcg$event, evmcmcg$param), mean)

pal <- gg_color_hue(2)
x <- 0:400 / 100

# pdf('../ev_est_gamma.pdf', width = 5, height = 3.5)
ggplot() +
  geom_histogram(aes(x = evdat$met_ndist, y = ..density..),
                 col = pal[1], fill = pal[1], alpha = 0.6) +
  geom_line(aes(x, dgamma(x, evpars['met a'], evpars['met b'])),
            col = pal[1], size = 1) +
  geom_histogram(aes(x = evdat$vert_ndist, y = ..density..),
                 col = pal[2], fill = pal[2], alpha = 0.6) +
  geom_line(aes(x, dgamma(x, evpars['vert a'], evpars['vert b'])),
            col = pal[2], size = 1) +
  xlab('Normalised distance') +
  ylab('Density') +
  xlim(0, 4)
# dev.off()

load('../outputs/sp2sp_mcmc_gathered.RData')

hupars <- by(humcmcg$value, paste(humcmcg$spto, humcmcg$param), mean)
yepars <- by(yemcmcg$value, paste(yemcmcg$spto, yemcmcg$param), mean)

remove(list = c('humcmcg', 'yemcmcg'))

load('../../03_stats_sp2sp/data/seed2sp_dist.RData')

yespl <- yedat$sp_to[!duplicated(yedat$sp_to)]

yepal <- gg_color_hue(length(yespl))
names(yepal) <- yespl

# pdf('../outputs/yeseed2sp_est_gamma.pdf', width = 5, height = 3.65)
for (sp in yespl) {
  alp <- yepars[paste(sp, 'a')]
  bet <- yepars[paste(sp, 'b')]
  
  ydat <- yedat[which(yedat$sp_to == sp), 'ndist_A']
  xl <- quantile(ydat, 0.99)
  
  x <- seq(0, xl, length.out = 1000)
  
  opl <- ggplot() +
    geom_histogram(aes(x = ydat, y = ..density..), alpha = 0.6,
                   col = yepal[sp], fill = yepal[sp]) +
    geom_line(aes(x, dgamma(x, alp, bet)), col = yepal[sp], size = 1) +
    xlim(0, xl) +
    labs(title = paste('Yeast phylome', sp)) +
    xlab('Normalised distance') +
    ylab('Density')
  
  print(opl)
}
# dev.off()


huspl <- hudat$sp_to[!duplicated(hudat$sp_to)]

hupal <- gg_color_hue(length(huspl))
names(hupal) <- huspl

# pdf('../outputs/huseed2sp_est_gamma.pdf', width = 5, height = 3.65)
for (sp in huspl) {
  alp <- hupars[paste(sp, 'a')]
  bet <- hupars[paste(sp, 'b')]
  
  ydat <- hudat[which(hudat$sp_to == sp), 'ndist_A']
  xl <- quantile(ydat, 0.99)
  
  x <- seq(0, xl, length.out = 1000)
  
  opl <- ggplot() +
    geom_histogram(aes(x = ydat, y = ..density..), alpha = 0.6,
                   col = hupal[sp], fill = hupal[sp]) +
    geom_line(aes(x, dgamma(x, alp, bet)), col = hupal[sp], size = 1) +
    xlim(0, xl) +
    labs(title = paste('Human phylome', sp)) +
    xlab('Normalised distance') +
    ylab('Density')
  
  print(opl)
}
# dev.off()