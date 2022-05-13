library(R2WinBUGS)

N <- 1000
nchains <- 3
niter <- 10000
thin <- 3
unifmax <- 100
seed <- 03820

model <- function() {
  # Likelihood
  for (i in 1:n) {
    y[i] ~ dgamma(a, b)
  }
  
  m <- a / b
  v <- a / pow(b, 2)
  mo <- (a - 1) / b
  
  # Priors
  a ~ dunif(0, unifmax)
  b ~ dunif(0, unifmax)
}

set.seed(seed)
y <- rgamma(N, 2, 4)

datlist <- list(y = y, n = N, unifmax = unifmax)

inits <- function() {
  list(a = runif(1, 0, 6), b = runif(1, 0, 6))
}

params <- c('a', 'b', 'm', 'v', 'mo')

post_raw <- bugs(model = model,
                 data = datlist,
                 inits = inits,
                 param = params,
                 n.chains = nchains,
                 n.burnin = round(niter * 0.1, 0),
                 n.thin = thin,
                 n.iter = niter,
                 WINE = '/usr/bin/wine',
                 bugs.directory = '/home/mgil/software/winbugs14_full_patched/WinBUGS14/',
                 debug = FALSE)

post <- mcmc.list()
for (i in 1:nchains) {
  post[[i]] <- mcmc(post_raw$sims.array[, i, ])
}

save(N, nchains, niter, thin, unifmax, y, post, post_raw,
     file = sprintf('../outputs_up/wbugs_it%sK_th%s_n%s_un%s_s%s.RData',
                    niter/1000, thin, N, unifmax, seed))
