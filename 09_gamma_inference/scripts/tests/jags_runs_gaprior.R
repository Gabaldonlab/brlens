library(rjags)

N <- 1000
nchains <- 3
niter <- 50000
thin <- 3
unifmax <- 100
seed <- 1

model_text <- 'model{
  # Likelihood
  for (i in 1:n) {
    y[i] ~ dgamma(a, b)
  }
  
  m <- a / b
  v <- a / b^2
  mo <- (a - 1) / b
  
  # Prior distributions
  a ~ dgamma(0.1, 0.1)
  b ~ dgamma(0.1, 0.1)
}'

set.seed(seed)
y <- rgamma(N, 2, 4)

datlist <- list(y = y, n = length(y))

model_jags <- jags.model(textConnection(model_text), 
                         data = datlist,
                         n.chains = nchains,
                         n.adapt = round(niter * 0.1, digits = 0))

post <- coda.samples(model_jags, 
                     variable.names = c('a', 'b', 'm', 'v', 'mo'), 
                     n.iter = niter,
                     thin = thin)

save(N, nchains, niter, thin, unifmax, y, post,
     file = sprintf('../outputs_up/jags_gaprior_it%sK_th%s_n%s_un%s_s%s.RData',
                    niter/1000, thin, N, unifmax, seed))
