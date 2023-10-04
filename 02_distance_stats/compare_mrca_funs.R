library(ggplot2)
library(tidyr)
library(ggpubr)
library(stringr)
library(see)

theme_set(theme_bw())

# Functions definition ----
tr_dat <- function(dat, nlabs) {
  ddat <- gather(dat[, c(1, grep('_dist', names(dat)))],
                 key = 'node', value = 'dist', -seed)
  ndat <- gather(dat[, c(1, grep('_ndist', names(dat)))],
                 key = 'node', value = 'ndist', -seed)
  
  ddat$node <- str_split(ddat$node, pattern = '_', simplify = TRUE)[, 1]
  ndat$node <- str_split(ndat$node, pattern = '_', simplify = TRUE)[, 1]
  
  mdat <- merge(ddat, ndat, by = c('seed', 'node'))
  mdat <- na.omit(mdat)
  
  mdat$node <- factor(mdat$node, levels = nlabs$Node, labels = nlabs$Abr)
  
  return(mdat)
}

# Loading data ----
new <- read.csv('../../01_distances/outputs_new/mammal_dist_new.csv', check.names = FALSE)
old <- read.csv('../../01_distances/outputs_old/mammal_dist_old.csv', check.names = FALSE)

nlabs <- read.csv('../data/node_labs.tsv', sep = '\t')
nlabs <- nlabs[order(nlabs$Node), ]
nlabs <- nlabs[-5, ]

gnew <- tr_dat(new, nlabs)
gold <- tr_dat(old, nlabs)

new_trees <- names(table(gnew$seed))
old_trees <- names(table(gold$seed))

# Filtering data ----
fgnew <- gnew[gnew$ndist < quantile(gnew$ndist, 0.99), ]
fgold <- gold[gold$ndist < quantile(gold$ndist, 0.99), ]

fnew_trees <- names(table(fgnew$seed))
fold_trees <- names(table(fgold$seed))

# save(file = '../data/dist_data.RData', fgnew, fgold)

# Plotting number of trees ----
# Raw data
ntrees <- data.frame(fun = c('New', 'Old'),
                     ntrees = c(length(new_trees), length(old_trees)),
                     check.names = FALSE)

a <- ggplot(ntrees, aes(reorder(fun, ntrees), ntrees, fill = fun)) +
  geom_bar(stat = 'identity') +
  xlab('MRCA function') +
  ylab('Number of trees') +
  geom_text(aes(y = ntrees / 2, x = fun, label = ntrees)) +
  scale_fill_okabeito()

tno <- merge(as.data.frame(table(gold$node)),
             as.data.frame(table(gnew$node)),
             by = 'Var1')
names(tno) <- c('Node', 'Old', 'New')

tno <- gather(tno, key = 'Method', value = 'Freq', -Node)

b <- ggplot(tno, aes(x = reorder(Node, -Freq), y = Freq, fill = Method)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  xlab('Node') +
  ylab('Number of trees') +
  geom_text(aes(y = Freq / 2, label = Freq),
            position = position_dodge(width = .9)) +
  scale_fill_okabeito()

# Filtered data
ntrees <- data.frame(fun = c('New', 'Old'),
                     ntrees = c(length(fnew_trees), length(fold_trees)),
                     check.names = FALSE)

c <- ggplot(ntrees, aes(reorder(fun, ntrees), ntrees, fill = fun)) +
  geom_bar(stat = 'identity') +
  xlab('MRCA function') +
  ylab('Number of trees') +
  geom_text(aes(y = ntrees / 2, x = fun, label = ntrees)) +
  scale_fill_okabeito()

tno <- merge(as.data.frame(table(fgold$node)),
             as.data.frame(table(fgnew$node)),
             by = 'Var1')
names(tno) <- c('Node', 'Old', 'New')

tno <- gather(tno, key = 'Method', value = 'Freq', -Node)

d <- ggplot(tno, aes(x = reorder(Node, -Freq), y = Freq, fill = Method)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  xlab('Node') +
  ylab('Number of trees') +
  geom_text(aes(y = Freq / 2, label = Freq),
            position = position_dodge(width = .9)) +
  scale_fill_okabeito()

ggarrange(a, b, c, d, common.legend = TRUE, widths = c(1, 2, 1, 2),
          labels = c('a)', 'b)', 'c)', 'd)'))

# Distribution plots ----
# Raw data, boxplots
a <- ggplot(gold, aes(y = dist, x = node, fill = node, colour = node)) +
  geom_boxplot(alpha = 0.6) +
  ylim(0, 2.5) +
  ylab('Raw distance') +
  xlab('Node')

b <- ggplot(gold, aes(y = ndist, x = node, fill = node, colour = node)) +
  geom_boxplot(alpha = 0.6) +
  ylim(0, 15) +
  ylab('Normalised distance') +
  xlab('Node')

c <- ggplot(gnew, aes(y = dist, x = node, fill = node, colour = node)) +
  geom_boxplot(alpha = 0.6) +
  ylim(0, 1) +
  ylab('Raw distance') +
  xlab('Node')

d <- ggplot(gnew, aes(y = ndist, x = node, fill = node, colour = node)) +
  geom_boxplot(alpha = 0.6) +
  ylim(0, 10) +
  ylab('Normalised distance') +
  xlab('Node')

ggarrange(a, b, c, d, common.legend = TRUE, align = 'hv',
          labels = 'auto')

# Raw data, density plots
a <- ggplot(gold, aes(x = dist, fill = node, colour = node)) +
  geom_density(alpha = 0.6) +
  xlim(0, 2.5) +
  xlab('Raw distance') +
  ylab('Density')

b <- ggplot(gold, aes(x = ndist, fill = node, colour = node)) +
  geom_density(alpha = 0.6) +
  xlim(0, 15) +
  xlab('Normalised distance') +
  ylab('Density')

c <- ggplot(gnew, aes(x = dist, fill = node, colour = node)) +
  geom_density(alpha = 0.6) +
  xlim(0, 1) +
  xlab('Raw distance') +
  ylab('Density')

d <- ggplot(gnew, aes(x = ndist, fill = node, colour = node)) +
  geom_density(alpha = 0.6) +
  xlim(0, 10) +
  xlab('Normalised distance') +
  ylab('Density')

ggarrange(a, b, c, d, common.legend = TRUE, align = 'hv',
          labels = 'auto')

# filtered data, boxplots
a <- ggplot(fgold, aes(y = dist, x = node, fill = node, colour = node)) +
  geom_boxplot(alpha = 0.6) +
  ylim(0, 4) +
  ylab('Raw distance') +
  xlab('Node')

b <- ggplot(fgold, aes(y = ndist, x = node, fill = node, colour = node)) +
  geom_boxplot(alpha = 0.6) +
  # ylim(0, 15) +
  ylab('Normalised distance') +
  xlab('Node')

c <- ggplot(fgnew, aes(y = dist, x = node, fill = node, colour = node)) +
  geom_boxplot(alpha = 0.6) +
  ylim(0, 1) +
  ylab('Raw distance') +
  xlab('Node')

d <- ggplot(fgnew, aes(y = ndist, x = node, fill = node, colour = node)) +
  geom_boxplot(alpha = 0.6) +
  # ylim(0, 10) +
  ylab('Normalised distance') +
  xlab('Node')

ggarrange(a, b, c, d, common.legend = TRUE, align = 'hv',
          labels = c('O', NA, 'N', NA))

# Filtered data, density plots
a <- ggplot(fgold, aes(x = dist, fill = node, colour = node)) +
  geom_density(alpha = 0.6) +
  # xlim(0, 2.5) +
  xlab('Raw distance') +
  ylab('Density')

b <- ggplot(fgold, aes(x = ndist, fill = node, colour = node)) +
  geom_density(alpha = 0.6) +
  # xlim(0, 15) +
  xlab('Normalised distance') +
  ylab('Density')

c <- ggplot(fgnew, aes(x = dist, fill = node, colour = node)) +
  geom_density(alpha = 0.6) +
  # xlim(0, 1) +
  xlab('Raw distance') +
  ylab('Density')

d <- ggplot(fgnew, aes(x = ndist, fill = node, colour = node)) +
  geom_density(alpha = 0.6) +
  # xlim(0, 10) +
  xlab('Normalised distance') +
  ylab('Density')

ggarrange(a, b, c, d, common.legend = TRUE, align = 'hv',
          labels = c('O', NA, 'N', NA))

library(ggridges)
ggplot(fgnew, aes(x = ndist, y = reorder(node, ndist),
                  fill = node, colour = node)) +
  geom_density_ridges(alpha = 0.6, show.legend = FALSE, scale = 5) +
  xlab('Normalised distance') +
  ylab('Density') +
  theme_ridges(center_axis_labels = TRUE)

