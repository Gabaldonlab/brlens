# Random phylome generation
# Moisès Bernabeu
# Barcelona, February 2022

# Libraries ----
library(ggtree)
library(treeio)
library(ape)
library(ggplot2)
library(gridExtra)
library(ggpubr)

theme_set(theme_bw())

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
    # rate <- rgamma(1, alpha, beta)
    rate <- sample(c(rgamma(1, 23, 12), rgamma(1, 1, 3)), 1)
    br_rates <- abs(rnorm(reftree$Nnode + length(reftree$tip.label) - 1,
                          rate, 0.1))
    nt$edge.length <- br_rates * reftree$edge.length
    nt <- list(nt)
    class(nt) <- 'multiPhylo'
    otrees <- c(otrees, nt)
  }
  return(otrees)
}

par(mfrow = c(1, 2))
x <- 0:400/100
plot(x, dgamma(x, 24, 10), type = 'l')
lines(x, dgamma(x, 2, 2))

plot(x, dgamma(x, 23, 12) + dgamma(x, 1, 3), type = 'l', ylab = 'density',
     xlab = 'rate')

rp <- ggplot() +
  geom_line(aes(x, dgamma(x, 23, 12) + dgamma(x, 1, 3))) +
  ylab('density') +
  xlab('rate')

# Reference tree ----
t_text <- '(((A:0.7,B:0.6):1.3,(C:0.3,D:0.5):2):1.2,E:3);'
t <- read.tree(text = t_text)

t_labs <- c('tau[A]', 'tau[B]', 'tau[C]', 'tau[D]', 'tau[E]', 'R',
            'tau[1]', 'tau[2]', 'tau[3]')
l_labs <- c('l[A]', 'l[B]', 'l[C]', 'l[D]', 'l[E]', 'R',
            'l[1]', 'l[2]', 'l[3]')

timetree <- ggtree(t) +
  geom_tiplab() +
  geom_text(aes(x = branch, label = t_labs),
            parse = TRUE, nudge_y = 0.15, size = 4)
timetree

ggtree(t) +
  geom_tiplab() +
  geom_label2(aes(x = branch, label = l_labs), parse = TRUE) +
  labs(title = 'Gene')

rp <- rp +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  ylab('Density') +
  xlab('Distance')

np <- ggplot() +
  geom_line(aes(0:400/100, dgamma(0:400/100, 3, 2))) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  ylab('Density') +
  xlab('Normalised distance')

# pdf('../outputs/densities_example.pdf', width = 8.3, height = 3)
ggarrange(rp, np, align = 'hv', labels = 'auto')
# dev.off()

rate <- ggplot() +
  geom_line(aes(0:400/100, y = dmixgamma(0:400/100))) +
  xlab('Rate') +
  ylab('Density')

# pdf('../compositions/p11_timetree.pdf', width = 7.5, height = 3.5)
ggarrange(timetree, rate, labels = 'auto', widths = c(1, 2))
# dev.off()

dmixgamma <- function(x, alpha = c(23, 1), beta = c(12, 3)) {
  k <- length(alpha)
  n <- length(x)
  rowSums(vapply(1:k, function(i) 1 / k * dgamma(x, alpha[i], beta[i]), numeric(n)))
}

# random generation
rmixgamma <- function(n, pi, alpha, beta) {
  k <- sample.int(length(pi), n, replace = TRUE, prob = pi)
  rgamma(n, alpha[k], beta[k])
}


set.seed(0233201384)
t <- rtree(20)

t$edge.length[which(t$edge.length < quantile(t$edge.length, 0.99))] <- t$edge.length[which(t$edge.length < quantile(t$edge.length, 0.99))] * 0.05

# pdf('../outputs/hr_tree.pdf', width = 2, height = 6)
ggtree(t)
# dev.off()

set.seed(0233201384)
t <- rtree(20)

t$edge.length[which(t$edge.length > quantile(t$edge.length, 0.05))] <- t$edge.length[which(t$edge.length > quantile(t$edge.length, 0.05))] * 5

# pdf('../outputs/lr_tree.pdf', width = 2, height = 6)
ggtree(t)
# dev.off()

p <- c()
for (i in 1:9) {
  p[[i]] <- ggtree(rtrees[[i]]) +
    geom_tiplab() +
    geom_label2(aes(x = branch, label = l_labs), parse = TRUE) +
    labs(title = paste0('GENE ', LETTERS[i], expression(paste0(' med(L[a])', LETTERS[i], '])'))), parse = TRUE)
}

do.call("grid.arrange", c(p, ncol=3))

# Random trees ----
rtrees <- rtreeset(t, n = 1000, mean_rate = 1)

# plot densitrees
ggdensitree(sample(rtrees, 100), alpha = 0.3, align.tips = TRUE) +
  geom_tiplab(align = TRUE)

# Write trees ----
write.tree(rtrees, file = '../outputs/rand_phylome.txt')

