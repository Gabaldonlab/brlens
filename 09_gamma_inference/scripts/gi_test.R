# Gamma inference - tests
# Moisès Bernabeu
# Barcelona, April 2022

library(ggplot2)
library(ggpubr)
library(rjags)
# library(DBDA2E)
library(coda)
library(ggmcmc)
library(tidyr)
library(wesanderson)
library(dplyr)
library(cumstats)

theme_set(theme_bw())

# Definitions ----
# Basic custom palette
pal <- c('steelblue', 'darkorange3', 'chartreuse3')

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


# Posterior gibbs ----
# With jags

model_text <- 'model{
  # Likelihood
  for (i in 1:n) {
    y[i] ~ dgamma(a, b)
  }
  
  m <- a / b
  v <- a / b^2
  mo <- (a - 1) / b
  
  # Prior distributions
  a ~ dunif(0, 100)
  b ~ dunif(0, 100)
}'

y <- rgamma(1000, 2, 4)

datlist <- list(y = y, n = length(y))

model_jags <- jags.model(textConnection(model_text), 
                         data = datlist,
                         n.chains = 3)

post <- coda.samples(model_jags, 
                     variable.names = c('a', 'b', 'm', 'v', 'mo'), 
                     n.iter = 1000)

summary(post)

post_list <- mcmc.list(post)
post_df <- as.data.frame(ggs(post_list))
post_df$Chain <- as.factor(post_df$Chain)

acarr <- autocorr(post, lags = 1:5)

autocorr.plot(post)

ggs_diagnostics(post_df)

ggplot(post_df, aes(x = value, colour = Chain)) +
  geom_density() +
  # geom_point(aes(y = 0), alpha = 0.1) +
  facet_wrap(~Parameter, scales = 'free') +
  scale_color_manual(values = wes_palette('FantasticFox1'))

ggplot(post_df, aes(x = Iteration, y = value, colour = Chain)) +
  geom_line() +
  facet_wrap(~Parameter, scales = 'free') +
  scale_color_manual(values = wes_palette('FantasticFox1'))

autocor_df <- data.frame()
for (chain in as.numeric(levels(post_df$Chain))) {
  for (param  in levels(post_df$Parameter)) {
    autocor_df <- rbind(autocor_df,
                        data.frame(Chain = paste(chain), Parameter = param,
                                   ac(post_df[which(post_df$Chain == chain &
                                                      post_df$Parameter == param),
                                              'value'],
                             nLags = 1000), stringsAsFactors = TRUE))
  }
}

ggplot(autocor_df, aes(x = Lag, y = Autocorrelation, colour = Chain)) +
  geom_line() +
  facet_grid(Parameter ~ Chain) +
  # xlim(0, 50) +
  scale_color_manual(values = wes_palette('FantasticFox1'))

# Calculate the mean of the chain
D <- post_df
dm.m <- D %>%
  dplyr::group_by(Parameter, Chain) %>%
  dplyr::summarize(m=mean(value))
# Calculate the running mean
# Force the object to be sorted by Parameter, and hence avoid 'rm' calculation
# to be wrong
dm.rm <- D %>%
  arrange(Parameter, Iteration) %>%
  group_by(Parameter, Chain) %>%
  mutate(rm = cumsum(value) / Iteration,
         sd = sqrt(cumvar(value)),
         up = rm + sd,
         do = rm - sd)

ggplot(dm.rm, aes(x = Iteration, y = rm, colour = as.factor(Chain))) +
  geom_line() +
  geom_line(aes(y = up), lty = 3, show.legend = FALSE) +
  geom_line(aes(y = do), lty = 3, show.legend = FALSE) +
  facet_grid(Parameter ~ Chain, scales = 'free') +
  geom_hline(aes(yintercept = m), dm.m, lty = 4, size = 0.3) +
  ylab('Running mean') +
  labs(colour = 'Chain')

ggs_running(post_df)
