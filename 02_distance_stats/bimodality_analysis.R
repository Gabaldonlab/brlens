library(tidyr)
library(ggfortify)
library(ggplot2)
theme_set(theme_bw())

load('../data/sp2sp_dat.RData')

plot(density(datc[which(datc$sp_to == 'PAPAN'), 'ndist']))
lines(density(datc[which(datc$sp_to == 'MACMU'), 'ndist']))
abline(v = 1.3)

# Selecting trees ----
plot(density(datc[which(datc$sp_to == 'PAPAN' & datc$ndist <= 1.5), 'ndist']))
lines(density(datc[which(datc$sp_to == 'MACMU' & datc$ndist <= 1.5), 'ndist']))

shrtpap <- datc[which(datc$sp_to == 'PAPAN' & datc$ndist <= 1.5), 'tree']
shrtpap <- shrtpap[!duplicated(shrtpap)]

shrtmac <- datc[which(datc$sp_to == 'MACMU' & datc$ndist <= 1.5), 'tree']
shrtmac <- shrtmac[!duplicated(shrtmac)]

sum(shrtmac %in% shrtpap) / length(shrtpap)
# 81.32% of the trees are lower for both species, it does not depend on the species

shrtall <- c(shrtmac, shrtpap)
shrtall <- shrtall[!duplicated(shrtall)]

plot(density(datc[which(datc$sp_to == 'PAPAN' & datc$ndist > 1.5), 'ndist']))
lines(density(datc[which(datc$sp_to == 'MACMU' & datc$ndist > 1.3), 'ndist']))

lngrpap <- datc[which(datc$sp_to == 'PAPAN' & datc$ndist > 1.5), 'tree']
lngrpap <- lngrpap[!duplicated(lngrpap)]

lngrmac <- datc[which(datc$sp_to == 'MACMU' & datc$ndist > 1.5), 'tree']
lngrmac <- lngrmac[!duplicated(lngrmac)]

sum(lngrpap %in% lngrmac) / length(lngrmac)
# 69.14% of the trees are lower for both species, it does not depend on the species

lngrall <- c(lngrmac, lngrpap)
lngrall <- lngrall[!duplicated(lngrall)]

# Generating the dataframe with all the distances and data ----
dat <- datc[which(datc$sp_to %in% c('PAPAN', 'MACMU')), c(1, 11, 10)]

treedists <- c()
for (tree in dat$tree[!duplicated(dat$tree)]) {
  ndf <- dat[which(dat$tree == tree), ]
  odf <- data.frame(tree = tree, MACMU = NA, PAPAN = NA)
  for (line in 1:dim(ndf)[1]) {
    odf[, ndf[line, 'sp_to']] <- ndf[line, 'ndist']
  }
  treedists <- rbind(treedists, odf)
}

# Relationships with tree variables ----
tstats <- read.csv('../../04_calc_dists/data/trees_stats.csv')
shrttdat <- tstats[which(tstats$tree %in% shrtall), ]
lngrtdat <- tstats[which(tstats$tree %in% lngrall), ]

stranged <- rbind(data.frame(shrttdat, group = 'shorter'),
                  data.frame(lngrtdat, group = 'longer'))

ggplot(stranged, aes(norm_width / tree_width, colour = group)) +
  geom_density()

ggplot(stranged, aes(norm_width / tree_width, colour = group)) +
  geom_density()

ggplot(stranged, aes(tree_skew / norm_skew, colour = group)) +
  geom_density() +
  xlim(-10, 10)

ggplot(stranged, aes(norm_sum / tree_sum, colour = group)) +
  geom_density() +
  labs(colour = 'Peak')

ggplot(stranged, aes(norm_leafno / tree_leafno, colour = group)) +
  geom_density() +
  labs(colour = 'Peak')

ndf <- merge(tstats, treedists, by = 'tree')
row.names(ndf) <- ndf$tree
ndf <- na.omit(ndf)
ndf$MACMU_lower <- ndf$MACMU <= 1.3
ndf$PAPAN_lower <- ndf$PAPAN <= 1.5
ndf[which(ndf$PAPAN_lower == ndf$MACMU_lower & ndf$MACMU_lower == TRUE), 'Group'] = 'Lower'
ndf[which(ndf$PAPAN_lower == ndf$MACMU_lower & ndf$MACMU_lower == FALSE), 'Group'] = 'Upper'
ndf[which(ndf$PAPAN_lower != ndf$MACMU_lower), 'Group'] = 'Disagree'

ggplot(ndf, aes(norm_leafno / norm_spno, MACMU, colour = MACMU_lower)) +
  geom_point()

ggplot(ndf, aes(norm_sum / tree_sum, PAPAN, colour = PAPAN_lower)) +
  geom_point()

pcadf <- na.omit(ndf[, -c(1, 20, 21, 22)])
plot(pcadf)

pca1 <- princomp(pcadf, cor = TRUE)

plot(pca1)

autoplot(pca1, data = ndf, colour = 'Group', loadings = TRUE, loadings.label = TRUE)

pca2df <- pcadf[, grep('tree', names(pcadf), value = TRUE)] / pcadf[, grep('norm', names(pcadf), value = TRUE)]
pca2df$MACMU <- pcadf$MACMU
pca2df$PAPAN <- pcadf$PAPAN

toget <- !is.infinite(pca2df$tree_skew)

pca2df <- pca2df[toget, ]

pca2 <- princomp(pca2df, cor = TRUE)

plot(pca2)

autoplot(pca2, data = ndf[toget, ], colour = 'Group', loadings = TRUE, loadings.label = TRUE)

pca3df <- ndf[-which(ndf$Group == 'Disagree'), ]
pca3df <- pca3df[, grep('tree_', names(pca3df), value = TRUE)] / pca3df[, grep('norm_', names(pca3df), value = TRUE)]
pca3df$MACMU <- ndf[-which(ndf$Group == 'Disagree'), ]$MACMU
pca3df$PAPAN <- ndf[-which(ndf$Group == 'Disagree'), ]$PAPAN

toget <- !is.infinite(pca3df$tree_skew)
pca3df <- pca3df[toget, ]

pca3 <- princomp(pca3df, cor = TRUE)

plot(pca3)

autoplot(pca3, data = ndf[-which(ndf$Group == 'Disagree'), ][toget, ],
         colour = 'Group', loadings = TRUE, loadings.label = TRUE)


torm <- c(ndf$tree[which(ndf$Group == 'Disagree')],
          row.names(pca3$scores)[which(pca3$scores[, 1] > 0.1)])

pca4df <- ndf[!(ndf$tree %in% torm), ]
pca4df <- pca4df[, grep('tree_', names(pca4df), value = TRUE)] / pca4df[, grep('norm_', names(pca4df), value = TRUE)]
pca4df$MACMU <- ndf[-which(ndf$Group == 'Disagree'), ]$MACMU
pca4df$PAPAN <- ndf[-which(ndf$Group == 'Disagree'), ]$PAPAN

toget <- !is.infinite(pca4df$tree_skew)
pca4df <- pca4df[toget, ]

pca4 <- princomp(pca4df, cor = TRUE)
plot(pca4)

autoplot(pca4, data = ndf[!(ndf$tree %in% torm), ][toget, ],
         colour = 'Group', loadings = TRUE, loadings.label = TRUE, alpha = 0.4)
