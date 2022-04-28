# Example plots for species to species distances
# Moisès Bernabeu
# Barcelona, April 2022

library(treeio)
library(ggtree)
library(gridExtra)
library(ggplot2)

# Example trees
set.seed(12031963)

tree <- rtree(20, tip.label = LETTERS[1:20])

d1920 <- c(20, 37, 38, 39, 19)
d0914 <- c(14, 9, 30, 28, 27, 26, 22, 21, 31, 32, 35)
noselnodes <- noselnodes <- c(1:39)[-c(d1920, d0914)]

cols <- c(rep('d(A, N)', length(d1920)),
          rep('d(C, P)', length(d0914)),
          rep(NA, length(noselnodes)))

d <- data.frame(node = c(d1920, d0914, noselnodes), group = cols)

tplot <- ggtree(tree, size = 1, aes(col = group)) %<+% d +
  geom_tippoint(aes(subset = node %in% c(19, 20)), size = 2, col = 'darkorange3') +
  geom_tippoint(aes(subset = node %in% c(9, 14)), size = 2, col = 'steelblue') +
  geom_tiplab(offset = 0.04, show.legend = FALSE) +
  # geom_tiplab(offset = 0.5, aes(label = node)) +
  scale_color_manual(values = c('d(C, P)' = 'steelblue',
                                'd(A, N)' = 'darkorange3'),
                     na.value = 'black') +
  labs(color = 'Distances')

# pdf('../msct_plots/tree.pdf', width = 4, height = 5.5)
tplot
# dev.off()

# Distance matrix
y <- read.csv('../../02_get_distances/outputs/0005_sptree_dist.csv')
y <- y[, c(4, 6, 10)]

nsp <- sum(!duplicated(y$from_sp))
sps <- y$from_sp[!duplicated(y$from_sp)]

m <- matrix(NA, nrow = nsp, ncol = nsp)
for (i in 1:nsp) {
  for (j in i:nsp) {
    if (sps[i] != sps[j]) {
      m[i, j] <- y[which((y$from_sp == sps[i] & y$to_sp == sps[j]) |
                           (y$from_sp == sps[j] & y$to_sp == sps[i])), 3]
      m[j, i] <- y[which((y$from_sp == sps[i] & y$to_sp == sps[j]) |
                           (y$from_sp == sps[j] & y$to_sp == sps[i])), 3]
    }
  }
}

# pdf('../msct_plots/heatmap.pdf', width = 5, height = 5)
heatmap(m, Rowv = NA, Colv = NA,
        labRow = LETTERS[1:20],
        labCol = LETTERS[1:20])
# dev.off()

# Phylome example
tplots <- c()
for (i in 1:20) {
  n <- rpois(1, 20)
  tree <- rtree(n)
  tplots[[i]] <- ggtree(tree)
}

# pdf('../msct_plots/phylome_2.pdf', width = 4, height = 7)
do.call(grid.arrange, tplots)
# dev.off()
