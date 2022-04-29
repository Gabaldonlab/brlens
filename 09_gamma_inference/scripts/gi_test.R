# Gamma inference - tests
# Moisès Bernabeu
# Barcelona, April 2022

library(ggplot2)
library(ggpubr)

theme_set(theme_bw())

# Likelihood
y <- rgamma(100, 2, 4)

alph <- 100:300/100
bet <- 300:600/100

like <- data.frame()
for (i in alph) {
  for (j in bet) {
    like <- rbind(like, data.frame(alpha = i, beta = j,
                                   like = prod(dgamma(y, i, j))))
  }
}

s100 <-  ggplot(like, aes(alpha, beta, fill = like)) +
  geom_tile() +
  labs(title = 'n = 100')


y <- rgamma(1000, 2, 4)

alph <- 100:300/100
bet <- 300:600/100

like <- data.frame()
for (i in alph) {
  for (j in bet) {
    like <- rbind(like, data.frame(alpha = i, beta = j,
                                   like = prod(dgamma(y, i, j))))
  }
}

s1000 <- ggplot(like, aes(alpha, beta, fill = like)) +
  geom_tile() +
  labs(title = 'n = 1000')

ggarrange(s100, s1000)


# Posterior gibbs
