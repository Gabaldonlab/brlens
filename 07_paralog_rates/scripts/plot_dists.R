library(tidyr)
library(ggplot2)
theme_set(theme_bw())

dat <- read.csv('../outputs/0076_dist.csv')

cols <- c('Mean' = 'steelblue', 'Median' = 'darkorange3')
ggplot(dat, aes(norm_D / norm_S, norm_median)) +
  geom_point() +
  geom_hline(yintercept = mean(dat$norm_median), linetype = 'Mean', col = 'steelblue') +
  geom_hline(yintercept = median(dat$norm_median), linetype = 'Median', col = 'darkorange3') +
  labs(color = 'Legend') +
  scale_color_manual(values = colors) +
  scale_linetype_manual(name = "limit", values = c(2, 4), 
                        guide = guide_legend(override.aes = list(color = c('steelblue', 'darkorange3'))))

gdupl <- gather(dat[, c('max_median', 'min_median')])

ggplot(gdupl, aes(y = value, x = key, col = key, fill = key)) +
  geom_boxplot(alpha = 0.6) +
  ylim(0, quantile(gdupl$value, 0.9))

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

