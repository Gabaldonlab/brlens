# Random phylome generation
# Moisès Bernabeu
# Barcelona, February 2022

# Libraries ----
library(ggtree)
library(treeio)
library(ape)

# Definitions ----
# Tree set definitions
rtreeset <- function(reftree, n = 100, mean_rate = 1) {
  otrees <- c()
  for (i in 1:n) {
    nt <- reftree
    
    # First we generate the protein rate by a gamma distribution with mean 1,
    # this is because of we want conservative trees
    alpha <- 2
    beta <- alpha / mean_rate
    rate <- rgamma(1, alpha, beta)
    br_rates <- abs(rnorm(reftree$Nnode + length(reftree$tip.label) - 1,
                          rate, 0.1))
    nt$edge.length <- br_rates * reftree$edge.length
    nt <- list(nt)
    class(nt) <- 'multiPhylo'
    otrees <- c(otrees, nt)
  }
  return(otrees)
}

# Reference tree ----
t_text <- '(((A:0.7,B:0.4):1,(C:0.12,D:0.5):2):1.2,E:1);'
t <- read.tree(text = t_text)

t_labs <- c('tau[A]', 'tau[B]', 'tau[C]', 'tau[D]', 'tau[E]', 'R',
            'tau[1]', 'tau[2]', 'tau[3]')

ggtree(t) +
  geom_tiplab() +
  geom_label2(aes(x = branch, label = t_labs), parse = TRUE)

# Random trees ----
rtrees <- rtreeset(t, n = 1000, mean_rate = 1)

# plot densitrees
ggdensitree(rtrees, alpha = 0.3)

# Write trees ----
write.tree(rtrees, file = '../outputs/rand_phylome.txt',tree.names = TRUE)
