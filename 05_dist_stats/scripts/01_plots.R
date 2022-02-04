library(ggplot2)

dat <- read.csv('../data/all_sumstats.csv', stringsAsFactors = TRUE)
str(dat)

ggplot(dat, aes(og_dist.Mean, og_sp.Mean)) +
  geom_point() +
  xlim(0, 10) +
  xlab('Outgroup-leaf mean distance') +
  ylab('Outgroup-leaf mean speciation events')

ggplot(dat, aes(og_dist.Median, og_sp.Median)) +
  geom_point() +
  xlim(0, 10) +
  xlab('Outgroup-leaf median distance') +
  ylab('Outgroup-leaf median speciation events')

ggplot(dat, aes(seed_dist.Mean, seed_sp.Mean)) +
  geom_point() +
  xlim(0, 10) +
  xlab('Seed-leaf mean distance') +
  ylab('Seed-leaf mean speciation events') +
  geom_smooth(method = 'lm')

ggplot(dat, aes(seed_dist.Median, seed_sp.Median)) +
  geom_point() +
  xlim(0, 10) +
  xlab('Seed-leaf median distance') +
  ylab('Seed-leaf median speciation events')

lm1 <- lm(dat$og_sp.Mean ~ dat$og_dist.Mean)
summary(lm1)

par(mfrow = c(2, 2))
plot(lm1)
dev.off()
