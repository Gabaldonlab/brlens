library(ggpubr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(ggtree)
library(treeio)

theme_set(theme_bw())

dat <- read.csv('../outputs/0076_dist.csv')

dat <- dat[-which(dat$vert_dist == dat$met_dist), ]

dat$wnwidthratio <- dat$whole_width / dat$norm_width

hist(wnwidthratio, breaks = 400, xlim = c(0, 100))
abline(v = quantile(wnwidthratio, 0.95))

dat <- dat[-which(wnwidthratio >= quantile(wnwidthratio, 0.95)), ]

rdistp <- ggplot(dat) +
  geom_density(aes(x = vert_dist, col = 'vertebrates')) +
  geom_density(aes(x = met_dist, col = 'metazoan')) +
  xlab('Seed to event raw distance') +
  labs(colour = 'Event') +
  ylab('Density')

# pdf('../outputs/raw_dist.pdf', width = 4.6, height = 2.4)
rdistp
# dev.off()

q99 <- quantile(c(dat$met_ndist, dat$vert_ndist), 0.99)

plot(dat$met_ndist, dat$whole_width)
text(dat$met_ndist[which(dat$met_ndist > q99)],
     dat$whole_width[which(dat$met_ndist > q99)] + 2,
     labels = dat$seed[which(dat$met_ndist > q99)])

# write.csv(dat$seed[which(dat$met_ndist > 20)], row.names = FALSE, file = '~/Desktop/miau.txt')

par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))
plot(dat$vert_ndist, xlab = 'Tree', ylab = 'Vertebrate distance')
plot(dat$met_ndist, xlab = 'Tree', ylab = 'Metazoan distance')
dev.off()

a <- ggplot() +
  geom_point(aes(dat$whole_width / dat$norm_width, dat$met_ndist / dat$met_dist)) +
  xlab('Whole tree width / Norm. clade width') +
  ylab('Metazoan normalised distance') +
  xlim(0, 300)

b <- ggplot() +
  geom_point(aes(dat$whole_width / dat$norm_width, dat$vert_ndist / dat$vert_dist)) +
  xlab('Whole tree width / Norm. clade width') +
  ylab('Vertebrate normalised distance') +
  xlim(0, 300)

ggarrange(a, b, hjust = 'hv')

# dat <- dat[-which(dat$vert_ndist > q09 & dat$met_ndist >= q09), ]

ndistp <- ggplot(dat) +
  geom_density(aes(x = vert_ndist, col = 'vertebrates')) +
  geom_density(aes(x = met_ndist, col = 'metazoan')) +
  xlim(0, quantile(c(dat$met_ndist, dat$vert_ndist), 0.9)) +
  xlab('Seed to event normalised distance') +
  labs(colour = 'Event')

dat[which(dat$met_ndist >= 20 | dat$vert_ndist >= 20), 1]

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
  ylim(0, quantile(grdat$distance, 0.9))

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

# save(rbp, nbp, ndistp, grdat, gndat, file = '~/Documents/research/working/dat.Rdata')

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

ggplot(dat, aes(x = whole_leafno, y = nfactor)) +
  geom_point() +
  xlim(0, 20)

data.frame(names(dat))
plot(dat[which(dat$met_dist > quantile(dat$met_dist, 0.9)), c(9:12, 19:20)])

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

# To look quitely!
datl1 <- dat[which(dat$nfactor < 1), ]
d4 <- datl1[, c(9:12, 19:20, 26)]

ps <- c()
for (v in 2:(length(d3))) {
  p <- ggplot(data.frame(x = d4[, v], y = log(d4[, 1])), aes(x, y)) +
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

trees <- read.tree('../outputs/')

library(stringr)
str_split(names(trees), ':', simplify = TRUE)[, 2]

p <- list()
i = 1
for (tree in trees) {
  p[[i]] <- ggtree(tree)
  i = i + 1
}

do.call("grid.arrange", c(p, ncol=7))
do.call("grid.arrange", p)

ggplot(dat, aes(n_dupl / n_sp, nfactor)) +
  geom_point() +
  ylim(0, 25) +
  xlab('Duplication rate') +
  ylab('Vertebrate norm. dist.')

a <- ggplot(dat, aes(n_dupl / n_sp, vert_ndist)) +
  geom_point() +
  ylim(0, 25) +
  xlab('Duplication rate') +
  ylab('Vertebrate norm. dist.')

b <- ggplot(dat, aes(n_dupl / n_sp, met_ndist)) +
  geom_point() +
  ylim(0, 25) +
  xlab('Duplication rate') +
  ylab('Metazoan norm. dist.')

pdf('../outputs/event_duprate.pdf', width = 9, height = 3)
ggarrange(a, b, align = 'hv', labels = 'auto')
dev.off()
