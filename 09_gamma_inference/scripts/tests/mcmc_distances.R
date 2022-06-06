library(rjags)

# https://sourceforge.net/p/mcmc-jags/discussion/610036/thread/585b0e4c/

mcmcfun <- function(spto, dat, nchains = 3, niter = 1000,
                    thin = 3, unifmax = 100) {
  require(rjags)
  
  model_text <- sprintf('model{
    # Likelihood
    for (i in 1:n) {
      y[i] ~ dgamma(a, b)
    }
    
    m <- a / b
    v <- a / b^2
    mo <- (a - 1) / b
    
    # Prior distributions
    a ~ dunif(0, %d)
    b ~ dunif(0, %d)
  }', unifmax, unifmax)
  
  sps <- dat$sp_to[!duplicated(dat$sp_to)]
  
  y <- dat[which(dat$sp_to == spto & dat$ndist_A != 0), 'ndist_A']
  N <- length(y)
  
  datlist <- list(y = y, n = N)
  
  model_jags <- jags.model(textConnection(model_text), 
                           data = datlist,
                           n.chains = nchains,
                           n.adapt = round(niter * 0.1, digits = 0))
  
  post <- coda.samples(model_jags, 
                       variable.names = c('a', 'b', 'm', 'v', 'mo'), 
                       n.iter = niter,
                       thin = thin)
  
  assign(sprintf('%s_mcmc', spto), post)
  
  save(list = c('N', 'nchains', 'niter', 'thin', 'unifmax', sprintf('%s_mcmc', spto)),
       file = sprintf('../outputs_up/jags_%s.RData', sprintf('%s_mcmc', spto)))
}

load('../../03_stats_sp2sp/data/seed2sp_dist.Rdata')

sps <- yedat$sp_to[!duplicated(yedat$sp_to)]

length(sps)

system.time(
  mcmcfun(sps[1], yedat, nchains = 3, niter = 1000, thin = 3, unifmax = 100)
)
