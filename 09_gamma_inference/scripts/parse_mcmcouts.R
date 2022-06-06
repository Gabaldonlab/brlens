# Parse MCMC of all the posteriors
# Moisès Bernabeu
# Barcelona, May 2022


library(stringr)
library(tidyr)
library(coda)

# Define functions
gathermcmc <- function(mcmcobj) {
  a <- as.matrix(get(mcmcobj))
  b <- gather(as.data.frame(a), key = 'param', value = 'value')
  b$spto <- str_split(mcmcobj, '_', simplify = TRUE)[, 1]
  return(b)
}

# Read data
hufiles <- list.files('/data/moises/cluster/brlens/09_gamma_inference/sp2sp_outputs/',
                      pattern = 'jags_hudat.*RData', full.names = TRUE)

yefiles <- list.files('/data/moises/cluster/brlens/09_gamma_inference/sp2sp_outputs/',
                      pattern = 'jags_yedat.*RData', full.names = TRUE)

event <- list.files('/data/moises/cluster/brlens/09_gamma_inference/outputs/',
                    pattern = 'RData', full.names = TRUE)

for (hufile in hufiles) {
  load(hufile)
}

humcmcg <- data.frame()
for (mcmcobj in ls(pattern = '_mcmc')) {
  print(mcmcobj)
  humcmcg <- rbind(humcmcg, gathermcmc(mcmcobj))
}

remove(list = ls(pattern = '_mcmc'))

for (yefile in yefiles) {
  load(yefile)
}

yemcmcg <- data.frame()
for (mcmcobj in ls(pattern = '_mcmc')) {
  print(mcmcobj)
  yemcmcg <- rbind(yemcmcg, gathermcmc(mcmcobj))
}

save(humcmcg, yemcmcg, file = '../outputs/mcmc_gathered.RData')

remove(list = ls(pattern = '_mcmc'))

for (ev in event) {
  load(ev)
}

evmcmcg <- data.frame()
for (mcmcobj in ls(pattern = '_mcmc')) {
  print(mcmcobj)
  evmcmcg <- rbind(evmcmcg, gathermcmc(mcmcobj))
}

save(evmcmcg, file = '../outputs/event_mcmc_gathered.RData')
