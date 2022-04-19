library(tidyr)
library(ggpubr)
library(ggplot2)
theme_set(theme_bw())

dat <- read.csv('../outputs/0076_dist.csv')

cols <- c('Mean' = 'steelblue', 'Median' = 'darkorange3')
ggplot(dat, aes(norm_D / norm_S, norm_median)) +
  geom_point() +
  geom_hline(yintercept = mean(dat$norm_median), lty = 2) +
  geom_hline(yintercept = median(dat$norm_median), lty = 4) +
  labs(color = 'Legend')

# Raw medians
meds <- dat[, c('max_median', 'min_median')]
summary(mmeds)
apply(meds, 2, sd)

gndupl <- gather(nmeds)

a <- ggplot(gdupl, aes(y = value, x = key, col = key, fill = key)) +
  geom_boxplot(alpha = 0.6) +
  ylim(0, quantile(gdupl$value, 0.9))

# Normalised
nmeds <- dat[, c('max_median', 'min_median')] / dat[, c('norm_median')]
summary(nmeds)
apply(nmeds, 2, sd)

gndupl <- gather(nmeds)

b <- ggplot(gndupl, aes(y = value, x = key, col = key, fill = key)) +
  geom_boxplot(alpha = 0.6) +
  ylim(0, quantile(gdupl$value, 0.9))

ggarrange(a, b, align = 'hv', common.legend = TRUE)

gdupl <- gather(dat[, c('max_mean', 'min_mean')])

ggplot(gdupl, aes(y = value, x = key, col = key, fill = key)) +
  geom_boxplot(alpha = 0.6) +
  ylim(0, quantile(gdupl$value, 0.9))

gdupl <- gather(dat[, c('max_skew', 'min_skew')])
ggplot(gdupl, aes(y = value, x = key, col = key, fill = key)) +
  geom_boxplot(alpha = 0.6)

# No var equal
var.test(dat$min_median, dat$max_median)

t.test(dat$min_median, dat$max_median, var.equal = FALSE)

gmeds <- gather(dat[, c('whole_D', 'min_median', 'max_median')],
                key = 'subtree', value = 'median',
                -whole_D)

dat <- read.csv('../outputs/0005_dist.csv')

ggplot(dat, aes(norm_median, norm_D / norm_S)) +
  geom_point()

gdupl <- gather(dat[, c('max_median', 'min_median')])

ggplot(gdupl, aes(y = value, x = key, col = key, fill = key)) +
  geom_boxplot(alpha = 0.6) +
  ylim(0, quantile(gdupl$value, 0.9)) +
  geom_hline(yu)

gdupl <- gather(dat[, c('max_mean', 'min_mean')])

ggplot(gdupl, aes(y = value, x = key, col = key, fill = key)) +
  geom_boxplot(alpha = 0.6) +
  ylim(0, quantile(gdupl$value, 0.9))

gdupl <- gather(dat[, c('max_skew', 'min_skew')])
ggplot(gdupl, aes(y = value, x = key, col = key, fill = key)) +
  geom_boxplot(alpha = 0.6)

# No var equal
var.test(dat$min_median, dat$max_median)

t.test(dat$min_median, dat$max_median, var.equal = FALSE)

gmeds <- gather(dat[, c('whole_D', 'min_median', 'max_median')],
                key = 'subtree', value = 'median',
                -whole_D)

