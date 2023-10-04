# Plot seed to species distance
# Moisès Bernabeu
# Barcelona, October 2022

# Loading libraries ---
library(ggplot2)
library(tidyr)
library(ggpubr)
library(forcats)

theme_set(theme_bw())

# Defining functions ----
distmode <- function(x) {
  distr <- density(x[x < quantile(x, 0.95, na.rm = TRUE)], na.rm = TRUE)
  y <- distr$x[which.max(distr$y)]
  return(y)
}

# Loading data ----
load('../data/sp2sp_dat.RData')

datc <- datc[which(datc$ndist <= quantile(datc$ndist, 0.9)), ]

ggplot(datc[which(datc$sp_to %in% c('PAPAN', 'MACMU')), ],
       aes(ndist, dup_count / (dup_count + sp_count), colour = sp_to)) +
  geom_point() +
  xlim(0, 1.5)

ggplot(datc[which(datc$sp_to %in% c('PAPAN', 'MACMU')), ],
       aes(dist, colour = sp_to)) +
  geom_density()

a <- datc[which(datc$sp_to %in% c('PAPAN', 'MACMU') & datc$ndist < 1.5), ]

ggplot(a, aes(dup_count / (dup_count + sp_count))) +
  geom_density()

b <- datc[which(datc$sp_to %in% c('PAPAN', 'MACMU') & datc$ndist >= 1.5), ]

ggplot(b, aes(dup_count / (dup_count + sp_count))) +
  geom_density()

ggplot(datc, aes(ndist, colour = sp_to)) +
  geom_density()

a <- ggplot(datc, aes(dist, colour = sp_to)) +
  geom_density(show.legend = FALSE) +
  facet_wrap(~fct_reorder(sp_to, dist, distmode), scales = 'free_y') +
  xlab('Raw distance') +
  ylab('Density')

b <- ggplot(datc, aes(ndist, colour = sp_to)) +
  geom_density(show.legend = FALSE) +
  facet_wrap(~fct_reorder(sp_to, ndist, distmode), scales = 'free_y') +
  xlab('Normalised distance') +
  ylab('Density')

ggarrange(a, b, labels = c('a)', 'b)'))

by(datc$sp_to, datc$tree, table)

# Bimodality origin
library(ggtree)

trees <- read.tree('../../04_calc_dists/data/mammal_trees_rooted.nwk')

a <- datc[which(datc$sp_to == 'PAPAN'), ]
b <- a$tree[which(a$ndist < 1.5)]
c <- a$tree[which(a$ndist >= 1.5)]

bim1 <- b[!duplicated(b)]
bim2 <- c[!duplicated(c)]

cat(bim1, sep = '\n', file = '../outputs/bimlth1_5_papan_seed.txt')
cat(bim2, sep = '\n', file = '../outputs/bimgth1_5_papan_seed.txt')

# pdf('../outputs/bimlth05_trees_papan.pdf', width = 5.7, height = 7.7)
for (tree in bim1) {
  aux <- a[which(a$tree == tree), ]
  df <- data.frame(tip = c(aux[1, 2], aux[1, 4]), lab = 'dist')
  p <- ggtree(trees[[tree]]) %<+% df +
    geom_tiplab(aes(colour = lab), show.legend = FALSE) +
    scale_color_manual(values = c('dist' = 'blue'), na.value = 'black')
  print(p)
}
# dev.off()

# pdf('../outputs/bimgth05_trees_papan.pdf', width = 5.7, height = 7.7)
for (tree in bim2) {
  aux <- a[which(a$tree == tree), ]
  df <- data.frame(tip = c(aux[1, 2], aux[1, 4]), lab = 'dist')
  p <- ggtree(trees[[tree]]) %<+% df +
    geom_tiplab(aes(colour = lab), show.legend = FALSE) +
    scale_color_manual(values = c('dist' = 'blue'), na.value = 'black')
  print(p)
}
# dev.off()

a <- datc[which(datc$sp_to == 'MACMU'), ]
b <- a$tree[which(a$ndist < 1.5)]
c <- a$tree[which(a$ndist >= 1.5)]

bim1 <- b[!duplicated(b)]
bim2 <- c[!duplicated(c)]

cat(bim1, sep = '\n', file = '../outputs/bimlth1_5_macmu_seed.txt')
cat(bim2, sep = '\n', file = '../outputs/bimgth1_5_macmu_seed.txt')

# pdf('../outputs/bimlth05_trees_macmu.pdf', width = 5.7, height = 7.7)
for (tree in bim1) {
  aux <- a[which(a$tree == tree), ]
  df <- data.frame(tip = c(aux[1, 2], aux[1, 4]), lab = 'dist')
  p <- ggtree(trees[[tree]]) %<+% df +
    geom_tiplab(aes(colour = lab), show.legend = FALSE) +
    scale_color_manual(values = c('dist' = 'blue'), na.value = 'black')
  print(p)
}
# dev.off()

# pdf('../outputs/bimgth05_trees_macmu.pdf', width = 5.7, height = 7.7)
for (tree in bim2) {
  aux <- a[which(a$tree == tree), ]
  df <- data.frame(tip = c(aux[1, 2], aux[1, 4]), lab = 'dist')
  p <- ggtree(trees[[tree]]) %<+% df +
    geom_tiplab(aes(colour = lab), show.legend = FALSE) +
    scale_color_manual(values = c('dist' = 'blue'), na.value = 'black')
  print(p)
}
# dev.off()

# cat(bim1, sep = '\n', file = '../outputs/bimlth05.txt')
# cat(bim2, sep = '\n', file = '../outputs/bimgth05.txt')

dat <- read.csv('../data/mammal_dist_new_mrca.csv', row.names = 1)

dat <- rbind(data.frame(dat[bim1, ], set = '< 1.5'),
             data.frame(dat[bim2, ], set = '>= 1.5'))

a <- ggplot(dat, aes(tree_D / (tree_S + tree_D), colour = set)) +
  geom_density() +
  scale_color_manual(values = c('< 1.5' = 'steelblue', '>= 1.5' = 'darkorange3'))

b <- ggplot(dat, aes(tree_width, colour = set)) +
  geom_density() +
  scale_color_manual(values = c('< 1.5' = 'steelblue', '>= 1.5' = 'darkorange3'))

c <- ggplot(dat, aes(tree_leafno - norm_leafno, colour = set)) +
  geom_density() +
  scale_color_manual(values = c('< 1.5' = 'steelblue', '>= 1.5' = 'darkorange3'))
  

ggarrange(a, b, c, nrow = 1, common.legend = 1)
