library(ggpubr)
library(ggplot2)
library(tidyr)
library(gridExtra)

theme_set(theme_bw())

dat <- read.csv('../outputs/0076_dist.csv')

dat <- dat[-which(dat$vert_dist == dat$met_dist), ]

rdistp <- ggplot(dat) +
  geom_density(aes(x = vert_dist, col = 'vertebrates')) +
  geom_density(aes(x = met_dist, col = 'metazoan')) +
  xlab('Seed to event raw distance') +
  labs(colour = 'Event')

# pdf('../outputs/raw_dist.pdf', width = 6.5, height = 3.5)
rdistp
# dev.off()

q09 <- quantile(c(dat$met_ndist, dat$vert_ndist), 0.9)

dat <- dat[-which(dat$vert_ndist > q09 & dat$met_ndist >= q09), ]

ndistp <- ggplot(dat) +
  geom_density(aes(x = vert_ndist, col = 'vertebrates')) +
  geom_density(aes(x = met_ndist, col = 'metazoan')) +
  xlim(0, quantile(c(dat$met_ndist, dat$vert_ndist), 0.9)) +
  xlab('Seed to event normalised distance') +
  labs(colour = 'Event')


# pdf('../outputs/norm_dist.pdf', width = 6.5, height = 3.5)
ndistp
# dev.off()

# pdf('../outputs/rn_dist.pdf', width = 7.8, height = 3.2)
ggarrange(rdistp, ndistp, labels = 'auto', common.legend = TRUE)
# dev.off()

rdat <- dat[, c(1, 3, 4)]
names(rdat) <- c('seed', 'Vertebrate', 'Metazoan')
grdat <- gather(rdat, key = 'group', value = 'distance', -c(seed))

rbp <- ggplot(grdat, aes(y = distance, x = group, color = group, fill = group)) +
  geom_boxplot(alpha = 0.6) +
  geom_violin(fill = NA) +
  ylab('Raw distance') +
  ylim(0, q09)

rbp

ndat <- dat[, c(1, 6, 7)]
names(ndat) <- c('seed', 'Vertebrate', 'Metazoan')
gndat <- gather(ndat, key = 'group', value = 'distance', -c(seed))

nbp <- ggplot(gndat, aes(y = distance, x = group,  color = group, fill = group)) +
  geom_boxplot(alpha = 0.6) +
  geom_violin(fill = NA) +
  ylab('Normalised distance') +
  ylim(0, quantile(gndat$distance, 0.9))

nbp

save(rbp, nbp, ndistp, grdat, gndat, file = '~/Documents/research/working/dat.Rdata')

# pdf('../outputs/rn_boxplots.pdf', width = 9, height = 3.5)
ggarrange(rbp, nbp, labels = 'auto', align = 'h', common.legend = TRUE)
# dev.off()

ggplot(dat, aes(x = n_dupl / n_sp, y = nfactor)) +
  geom_point() +
  geom_hline(yintercept = median(dat$nfactor), col = 'steelblue', lty = 4) +
  geom_hline(yintercept = mean(dat$nfactor), col = 'darkorange3', lty = 4)

ggplot(dat, aes(x = whole_width, y = nfactor)) +
  geom_point() +
  xlim(0, 20)

ggplot(dat[which(dat$nfactor > 4.5), ], aes(x = whole_mean, y = nfactor)) +
  geom_point()

d2 <- dat[which(dat$nfactor > 4.5), ]

as.data.frame(names(d2))

plot(d2[, c(9:11, 20:23, 26, 27)], pch = 20)

data.frame(names(dat))
d3 <- dat[which(dat$nfactor > 4.5), c(9:12, 19:20, 26)]

ps <- c()
for (v in 2:(length(d3))) {
  p <- ggplot(data.frame(x = d3[, v], y = d3[, 1]), aes(x, y)) +
    geom_point() +
    ylab('Normalising factor') +
    xlab(names(d3)[v])
  ps[[v - 1]] <- p
}

do.call("grid.arrange", c(ps, ncol=3))

ks.test(dat$vert_ndist, y = 'pnorm')
ks.test(dat$met_ndist, y = 'pnorm')
shapiro.test(sample(dat$vert_ndist, 5000))
shapiro.test(sample(dat$met_ndist, 5000))
var.test(gndat$distance ~ gndat$group)
