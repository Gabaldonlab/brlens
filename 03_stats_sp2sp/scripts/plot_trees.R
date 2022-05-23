# Plotting species trees
# Moisès Bernabeu
# March, 2022

library(ggtree)
library(treeio)
library(ggpubr)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

tree <- read.tree('../../02_get_distances/data/0005_sptree.nwk')
tdat <- read.csv('../../02_get_distances/data/0005_norm_groups.csv')

yetree <- ggtree(tree) %<+% tdat +
  geom_tiplab(aes(label = Species.Name, colour = Normalising.group)) +
  geom_tree() +
  xlim(0, 3.5) +
  scale_color_manual(values = c('A' = 'steelblue'), na.value = 'black')
yetree

tree <- read.tree('../../02_get_distances/data/0076_sptree.nwk')
tdat <- read.csv('../../02_get_distances/data/0076_norm_groups.csv')

hutree <- ggtree(tree) %<+% tdat +
  geom_tiplab(aes(label = Species.Name, colour = Normalising.group)) +
  geom_tree() +
  xlim(0, 12) +
  scale_color_manual(values = c('A' = 'steelblue'), na.value = 'black') +
  geom_nodepoint(aes(subset = node %in% 51), size = 2, col = gg_color_hue(2)[2]) + # Vertebrates
  geom_nodepoint(aes(subset = node %in% 48), size = 2, col = gg_color_hue(1)) # Metazoans
hutree

pdf('../outputs/trees_plot.pdf', width = 9, height = 6.40)
ggarrange(yetree, hutree, common.legend = TRUE)
dev.off()

tree <- read.tree('../../02_get_distances/data/qfo78_sp_tree.txt')
tdat <- read.csv('../../02_get_distances/data/qfo_78_norm_groups.csv')
tdat$Proteome_ID <- tdat$Proteome

qfotree <- ggtree(tree) %<+% tdat +
  geom_tiplab(aes(colour = Normalising.group), size = 3) +
  geom_tree() +
  xlim(0, 20)
qfotree

