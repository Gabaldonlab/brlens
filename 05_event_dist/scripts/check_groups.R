library(treeio)
library(ggtree)

dat <- read.csv('../outputs/0076_dist.csv')
full_trees <- read.csv('../data/0076_108.txt', sep = '\t', header = FALSE)

tree_set <- which(dat$vert_dist > dat$met_dist)

dat <- dat[tree_set, ]
full_trees <- full_trees[tree_set, ]

'%!in%' <- function(x,y)!('%in%'(x,y))

for (i in 1:length(tree_set)) {
  t <- read.tree(text = full_trees[i, 4])
  tp <- ggtree(t) +
    geom_tiplab()
  
  write(file = 'tmp_tree.txt', dat[i, 'vert_node'])
  vt <- read.nhx('tmp_tree.txt')
  vtp <- ggtree(vt) +
    geom_tiplab(aes(col = Vertebrate))
  
  write(file = 'tmp_tree.txt', dat[i, 'met_node'])
  mt <- read.nhx('tmp_tree.txt')
  mtp <- ggtree(mt) +
    geom_tiplab(aes(col = Metazoan))
  
  df <- data.frame(taxa = t$tip.label)
  df[which(df$taxa %in% vt@phylo$tip.label & df$taxa %in% mt@phylo$tip.label), 'Group'] = 'vertebrate'
  df[which(df$taxa %in% mt@phylo$tip.label & df$taxa %!in% vt@phylo$tip.label), 'Group'] = 'metazoan'
  df[which(df$taxa == full_trees[i, 1]), 'Group'] <- 'SEED'
  
  ggtree(t) %<+% df +
    geom_tiplab(aes(color = Group))
}
