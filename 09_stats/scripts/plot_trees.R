# Plotting species trees
# Moisès Bernabeu
# March, 2022

library(ggtree)
library(treeio)
library(ggpubr)

tree <- read.tree('../../11_cladenorm/data/0005_sptree.nwk')
tdat <- read.csv('../../11_cladenorm/data/0005_norm_groups.csv')

yetree <- ggtree(tree) %<+% tdat +
  geom_tiplab(aes(colour = Normalising.group)) +
  geom_tree() +
  xlim(0, 2.5)
yetree

tree <- read.tree('../../11_cladenorm/data/0076_sptree.nwk')
tdat <- read.csv('../../11_cladenorm/data/0076_norm_groups.csv')

hutree <- ggtree(tree) %<+% tdat +
  geom_tiplab(aes(colour = Normalising.group)) +
  geom_tree() +
  xlim(0, 8)
hutree

ggarrange(yetree, hutree, common.legend = TRUE)

tree <- read.tree('../../11_cladenorm/data/qfo78_sp_tree.txt')
tdat <- read.csv('../../11_cladenorm/data/qfo_78_norm_groups.csv')
tdat$Proteome_ID <- tdat$Proteome

qfotree <- ggtree(tree) %<+% tdat +
  geom_tiplab(aes(colour = Normalising.group), size = 3) +
  geom_tree() +
  xlim(0, 20)
qfotree

