# Descriptive statistics and data filtering
# Moisès Bernabeu
# Sant Feliu de Guíxols, September 2022

# Loading libraries ----
library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)
library(modeest)
library(ggpubr)
library(see)

theme_set(theme_bw())

# Loading data ----
dat <- read.csv('../../01_distances/outputs_new/mammal_dist_new.csv', check.names = FALSE)
nlabs <- read.csv('../data/node_labs.tsv', sep = '\t')
nlabs <- nlabs[order(nlabs$Node), ]
nlabs <- nlabs[-5, ]

ddat <- gather(dat[, c(1, grep('_dist', names(dat)))],
               key = 'node', value = 'dist', -seed)
ndat <- gather(dat[, c(1, grep('_ndist', names(dat)))],
               key = 'node', value = 'ndist', -seed)

ddat$node <- str_split(ddat$node, pattern = '_', simplify = TRUE)[, 1]
ndat$node <- str_split(ndat$node, pattern = '_', simplify = TRUE)[, 1]

mdat <- merge(ddat, ndat, by = c('seed', 'node'))
mdat <- na.omit(mdat)

mdat$node <- factor(mdat$node, levels = nlabs$Node, labels = nlabs$Abr)

# mdat$seed[which(!duplicated(mdat$seed))]

# Basic density plots ----
ggplot(mdat, aes(dist, colour = node)) +
  geom_density()

ggplot(mdat, aes(ndist, colour = node)) +
  geom_density()

# Tree filtering ----
# Diagnostic plots
a <- ggplot(dat, aes(wdth_ratio)) +
  geom_density() +
  xlab('Width ratio') +
  ylab('Density')

b <- ggplot(dat, aes(tree_leafno, norm_leafno)) +
  geom_point() +
  xlab('Tree leaf no.') +
  ylab('Norm. clade. leaf no.')

c <- ggplot(dat, aes(wdth_ratio, norm_leafno / tree_leafno)) +
  geom_point() +
  xlab('Width ratio') +
  ylab('Tree and norm. clade leaf ratio')

# pdf('../outputs/raw_diagnostics.pdf', width = 12, height = 3, onefile = FALSE)
ggarrange(a, b, c, align = 'hv', nrow = 1, labels = 'auto')
# dev.off()

# Filtering trees
large_dists <- mdat[which(mdat$ndist > quantile(mdat$ndist, 0.99)), 'seed']
large_dists <- large_dists[!duplicated(large_dists)]

# cat(large_dists, sep = '\n', file = '../outputs/large_dist_trees_old_mrca.txt')

fdat <- dat[-which(dat$seed %in% large_dists), ]
fmdat <- mdat[-which(mdat$seed %in% large_dists), ]

# Diagnostic plots
a <- ggplot(fdat, aes(wdth_ratio)) +
  geom_density() +
  xlab('Width ratio') +
  ylab('Density')

b <- ggplot(fdat, aes(tree_leafno, norm_leafno)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  xlab('Tree leaf no.') +
  ylab('Norm. clade. leaf no.')

c <- ggplot(fdat, aes(wdth_ratio, norm_leafno / tree_leafno)) +
  geom_point() +
  xlab('Width ratio') +
  ylab('Tree and norm. clade leaf ratio')

# pdf('../outputs/filt_diagnostics.pdf', width = 12, height = 3, onefile = FALSE)
ggarrange(a, b, c, align = 'hv', nrow = 1, labels = 'auto')
# dev.off()

# Plotting filtered data ----
dbp <- ggplot(fmdat, aes(node, dist, colour = node)) +
  geom_boxplot(aes(fill = node), alpha = 0.6) +
  ylim(0, 1) +
  ylab('Distance') +
  xlab('Node')

ndbp <- ggplot(fmdat, aes(node, ndist, colour = node)) +
  geom_boxplot(aes(fill = node), alpha = 0.6) +
  # ylim(0, 25) +
  ylab('Normalised distance') +
  xlab('Node')

# pdf('../outputs/dists_boxplots.pdf', width = 7, height = 3, onefile = FALSE)
ggarrange(dbp, ndbp, align = 'hv', legend = 'bottom', common.legend = TRUE,
          labels = 'auto')
# dev.off()

treeno <- as.data.frame(table(fmdat$node))

ggplot(treeno, aes(reorder(Var1, -Freq), Freq)) +
  geom_bar(stat = 'identity', fill = 'steelblue') +
  xlab('Node') +
  ylab('Number of trees') +
  geom_text(aes(y = Freq / 2, label = Freq))

ddens <- ggplot(fmdat, aes(dist, colour = node)) +
  geom_density() +
  xlim(0, 2.5) +
  ylab('Density') +
  xlab('Distance') +
  scale_color_okabeito()

nddens <- ggplot(fmdat, aes(ndist, colour = node)) +
  geom_density() +
  xlim(0, 7.5) +
  ylab('Density') +
  xlab('Normalised distance') +
  scale_color_okabeito()

# pdf('../outputs/dists_density.pdf', width = 7, height = 3, onefile = FALSE)
ggarrange(ddens, nddens, align = 'hv', legend = 'bottom', common.legend = TRUE,
          labels = c('a)', 'b)'))
# dev.off()

# Median's tests Mann-Withney's U test ----
wilcox.test(fmdat[which(fmdat$node == 'T'), 'ndist'],
            fmdat[which(fmdat$node == 'Pl'), 'ndist'], paired = FALSE,
            alternative = 'greater')
wilcox.test(fmdat[which(fmdat$node == 'Pl'), 'ndist'],
            fmdat[which(fmdat$node == 'B'), 'ndist'], paired = FALSE,
            alternative = 'greater')
wilcox.test(fmdat[which(fmdat$node == 'B'), 'ndist'],
            fmdat[which(fmdat$node == 'P'), 'ndist'], paired = FALSE,
            alternative = 'greater')

