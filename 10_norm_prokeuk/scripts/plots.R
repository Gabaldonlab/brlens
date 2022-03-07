library(ggplot2)
library(ggpubr)

theme_set(theme_bw())

dat <- rbind(data.frame(read.csv('../outputs/0739_dist_noh.csv'),
                        'phylome' = 'human'),
             data.frame(read.csv('../outputs/0907_dist_noh.csv'),
                        'phylome' = 'yeast'))

ggplot(dat, aes(med / median(med), col = subtree)) +
  geom_density() +
  geom_vline(xintercept = 1, col = 'darkgrey', lty = 4) +
  xlim(0, 15) +
  facet_wrap(~phylome + norm_fact, scales = 'free')
