# MCMC sampling for data
# Moisès Bernabeu
# Barcelona, May 2022

# Loading sampling function
source('mcmc_sampling_function.R')

# Loading data
load('../../03_stats_sp2sp/data/seed2sp_dist.Rdata')

# Getting species list
sps <- yedat$sp_to[!duplicated(yedat$sp_to)]
length(sps)

# Executing sampling
system.time(
  mcmcout <- mcmcfun(sps[1], yedat, nchains = 3, niter = 100000, thin = 3, unifmax = 100)
)

# Plotting posterior distributions
pdf(sprintf('../outputs/jags_%s_plots.pdf', sps[10]), width = 6.3, height = 7)
plot(mcmcout)
dev.off()