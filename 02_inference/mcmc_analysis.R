# Functions to analyse MCMC
# Moisès Bernabeu
# Barcelona, Janyary 2022

# Requiring packages ----
require(ggmcmc)
require(tidyr)
require(dplyr)
require(cumstats)
require(coda)
require(ggpubr)

theme_set(theme_bw())

# Defining function ----
get_ac_df <- function(mcmc_df, nLags = 50) {
  # Retrieve autocorrelation per parameter and chain
  wc.ac <- mcmc_df %>%
    group_by(Parameter, Chain) %>%
    do(ac(.$value, nLags))
  
  return(wc.ac)
}


get_full_means <- function(mcmc_df) {
  # Estimates of the mean by iteration
  dm.m <- mcmc_df %>%
    group_by(Parameter, Chain) %>%
    summarize(m = mean(value))
  
  return(dm.m)
}

get_running <- function(mcmc_df) {
  # Calculate the running mean
  # Force the object to be sorted by Parameter, and hence avoid 'rm' calculation
  # to be wrong
  dm.rm <- mcmc_df %>%
    arrange(Parameter, Iteration) %>%
    group_by(Parameter, Chain) %>%
    mutate(rm = cumsum(value) / Iteration,
           sd = sqrt(cumvar(value)),
           up = rm + sd,
           do = rm - sd)
  
  return(dm.rm)
}

plot_diag <- function(y, distr, params, title) {
  # Plots the inferred distribution against the data histogram, the
  # quantile-quantile plot and the CDF for a given distribution parameters
  # samples returning a list of plots.
  dat <- data.frame(y = y)
  
  x <- seq(range(dat)[1], range(dat)[2],  0.001)
  
  if (distr == 'gamma') {
    ytheo = dgamma(x, params[1], params[2])
    lcol <- 'steelblue'
    linety <- 1
    qq <- geom_qq(distribution = qgamma, dparams = params, colour = lcol)
    li <- geom_line(aes(x = x, y = pgamma(x, params[1], params[2])),
                    colour = lcol, lty = linety, size = 0.75)
    
  } else if (distr == 'normal') {
    ytheo = dnorm(x, params[1], params[2])
    lcol <- 'darkorange3'
    linety <- 4
    qq <- geom_qq(distribution = qnorm, dparams = params, colour = lcol)
    li <- geom_line(aes(x = x, y = pnorm(x, params[1], params[2])),
                    colour = lcol, lty = linety, size = 0.75)
  }
  
  infdat <- data.frame(x = x, y = ytheo)
  
  a <- ggplot(dat) +
    geom_histogram(aes(x = y, y = after_stat(density)), alpha = 0.4, colour = 'black') +
    geom_line(data = infdat, aes(x, y), size = 0.75, col = lcol, lty = linety) +
    ylab('Density') +
    xlab('Distance') +
    labs(colour = 'Inferred distribution', lty = 'Inferred distribution')
  
  b <- ggplot(dat, aes(sample = y)) +
    geom_abline(slope = 1) +
    qq +
    xlab('Theoretical') +
    ylab('Observed')
  
  c <- ggplot() +
    stat_ecdf(data = dat, aes(y), geom = 'point') +
    li +
    ylab('CDF') +
    xlab('Distance')
  
  p <- ggarrange(a, b, c, common.legend = TRUE, nrow = 1, align = 'hv')
  p <- annotate_figure(p, top = text_grob(title, color = 'black',
                                          face = 'bold', size = 14))
  print(p)
}

plot_diag_2distr <- function(y, gmcmcout, nmcmcout, title) {
  # Plots the inferred distribution against the data histogram, the
  # quantile-quantile plot and the CDF for both Gamma and Normal distribution
  # samples returning a plot with the 3 subplots.
  dat <- data.frame(y = y)
  
  x <- seq(range(dat)[1], range(dat)[2],  0.001)
  
  dparams <- list('Gamma' = summary(gmcmcout)$statistics[1:2, 'Mean'],
                  'Normal' = summary(nmcmcout)$statistics[1:2, 'Mean'])
  names(dparams[['Gamma']]) <- NULL
  names(dparams[['Normal']]) <- NULL
  
  infdat <- data.frame(x = c(x, x),
                       y = c(dgamma(x, dparams[['Gamma']][1], dparams[['Gamma']][2]),
                             dnorm(x, dparams[['Normal']][1], dparams[['Normal']][2])),
                       distr = c(rep('Gamma', length(x)),
                                 rep('Normal', length(x))))
  
  a <- ggplot(dat) +
    geom_histogram(aes(x = y, y = after_stat(density)), alpha = 0.4, colour = 'black') +
    geom_line(data = infdat, aes(x, y, colour = distr, lty = distr), size = 0.75) +
    scale_colour_manual(values = c('Gamma' = 'steelblue', 'Normal' = 'darkorange3')) +
    scale_linetype_manual(values = c('Gamma' = 1, 'Normal' = 4)) +
    ylab('Density') +
    xlab('Distance') +
    labs(colour = 'Inferred distribution', lty = 'Inferred distribution')
  
  b <- ggplot(dat, aes(sample = y)) +
    geom_abline(slope = 1) +
    geom_qq(distribution = qnorm, dparams = dparams$Normal, colour = 'darkorange3') +
    geom_qq(distribution = qgamma, dparams = dparams$Gamma, colour = 'steelblue') +
    xlab('Theoretical') +
    ylab('Observed')
  
  c <- ggplot() +
    stat_ecdf(data = dat, aes(y), geom = 'point') +
    geom_line(aes(x = x, y = pnorm(x, dparams$Normal[1], dparams$Normal[2])),
              colour = 'darkorange3', lty = 4, size = 0.75) +
    geom_line(aes(x = x, y = pgamma(x, dparams$Gamma[1], dparams$Gamma[2])),
              colour = 'steelblue', lty = 1, size = 0.75) +
    ylab('CDF') +
    xlab('Distance')
  
  p <- ggarrange(a, b, c, common.legend = TRUE, nrow = 1, align = 'hv')
  p <- annotate_figure(p, top = text_grob(title, color = 'black',
                                          face = 'bold', size = 14))
  return(p)
}

plot_diag_2distr_list <- function(y, gmcmcout, nmcmcout, title) {
  # Plots the inferred distribution against the data histogram, the
  # quantile-quantile plot and the CDF for both Gamma and Normal distribution
  # samples returning a list of plots.
  dat <- data.frame(y = y)
  
  x <- seq(range(dat)[1], range(dat)[2],  0.001)
  
  dparams <- list('Gamma' = summary(gmcmcout)$statistics[1:2, 'Mean'],
                  'Normal' = summary(nmcmcout)$statistics[1:2, 'Mean'])
  names(dparams[['Gamma']]) <- NULL
  names(dparams[['Normal']]) <- NULL
  
  infdat <- data.frame(x = c(x, x),
                       y = c(dgamma(x, dparams[['Gamma']][1], dparams[['Gamma']][2]),
                             dnorm(x, dparams[['Normal']][1], dparams[['Normal']][2])),
                       distr = c(rep('Gamma', length(x)),
                                 rep('Normal', length(x))))
  
  a <- ggplot(dat) +
    geom_histogram(aes(x = y, y = after_stat(density)), alpha = 0.4, colour = 'black') +
    geom_line(data = infdat, aes(x, y, colour = distr, lty = distr), size = 0.75) +
    scale_colour_manual(values = c('Gamma' = 'steelblue', 'Normal' = 'darkorange3')) +
    scale_linetype_manual(values = c('Gamma' = 1, 'Normal' = 4)) +
    ylab('Density') +
    xlab('Distance') +
    labs(colour = 'Inferred distribution', lty = 'Inferred distribution', title = title)
  
  b <- ggplot(dat, aes(sample = y)) +
    geom_abline(slope = 1) +
    geom_qq(distribution = qnorm, dparams = dparams$Normal, colour = 'darkorange3') +
    geom_qq(distribution = qgamma, dparams = dparams$Gamma, colour = 'steelblue') +
    xlab('Theoretical') +
    ylab('Observed') +
    labs(title = title)
  
  c <- ggplot() +
    stat_ecdf(data = dat, aes(y), geom = 'point') +
    geom_line(aes(x = x, y = pnorm(x, dparams$Normal[1], dparams$Normal[2])),
              colour = 'darkorange3', lty = 4, size = 0.75) +
    geom_line(aes(x = x, y = pgamma(x, dparams$Gamma[1], dparams$Gamma[2])),
              colour = 'steelblue', lty = 1, size = 0.75) +
    ylab('CDF') +
    xlab('Distance') +
    labs(title = title)
  
  p <- list(a + theme(legend.position = 'none'), b, c)
  
  return(p)
}

mcmc_to_df <- function(mcmc) {
  # Converts an MCMC list object into a dataframe to be used in ggplot
  post_list <- mcmc.list(mcmc)
  post_df <- as.data.frame(ggs(post_list))
  post_df$Chain <- as.factor(post_df$Chain)
  
  return(post_df)
}

ggtraces <- function(mcmc_df, title = '') {
  # Plots the MCMC traces, checking whether the file is an MCMC, if it is, the
  # function will convert it using mcmc_to_df
  if (class(mcmc_df) %in% c('mcmc.list', 'mcmc')) {
    mcmc_df <- mcmc_to_df(mcmc_df)
  }
  
  gtraces <- ggplot(mcmc_df, aes(x = Iteration, y = value, colour = Chain)) +
    geom_line(size = 0.2, aes(lty = Chain)) +
    facet_wrap(~Parameter, scales = 'free', ncol = 1) +
    labs(title = title)
  
  return(gtraces)
}

ggmcmcdens <- function(mcmc_df, title = '') {
  # Plots the MCMC densities, checking whether the file is an MCMC, if it is, the
  # function will convert it using mcmc_to_df
  if (class(mcmc_df) %in% c('mcmc.list', 'mcmc')) {
    mcmc_df <- mcmc_to_df(mcmc_df)
  }

  ghists <- ggplot(mcmc_df, aes(x = value, colour = Chain)) +
    geom_density() +
    facet_wrap(~Parameter, scales = 'free', ncol = 1) +
    labs(title = title)
  
  return(ghists)
}

ggmcmcsummary <- function(mcmc_df, cl_mcmc_df, title = '') {
  # Plots the MCMC densities, checking whether the file is an MCMC, if it is, the
  # function will convert it using mcmc_to_df
  p <- ggarrange(ggtraces(mcmc_df), ggtraces(cl_mcmc_df), ggmcmcdens(cl_mcmc_df),
                 ncol = 3, align = 'hv', common.legend = TRUE,
                 legend = 'right')
  if (title != '') {
    p <- annotate_figure(p, top = text_grob(title, face = "bold", size = 14))
  }

  return(p)
}

numeric_diagnostics <- function(mcmc, mcmc_df = '') {
  if (mcmc_df == '') {
    mcmc_df <- mcmc_to_df(mcmc)
  }
  
  diagdf <- as.data.frame(ggs_diagnostics(mcmc_df))
  summ <- summary(mcmc)
  summdf <- cbind(summ$statistics[, 1:2],
                  'R' = diagdf[which(diagdf$Diagnostic == 'Rhat'), 4],
                  'Effective size' = effectiveSize(mcmc))
  
  return(summdf)
}

ggacplot <- function(mcmc_df, title = '') {
  if (class(mcmc_df) %in% c('mcmc.list', 'mcmc')) {
    mcmc_df <- mcmc_to_df(mcmc_df)
  }

  post_ac <- get_ac_df(mcmc_df, 100)
  ac <- ggplot(post_ac, aes(x = Lag, y = Autocorrelation, colour = Chain)) +
    geom_line() +
    facet_grid(Parameter ~ Chain) +
    labs(title = title)
  
  return(ac)
}

graphic_diagnostics <- function(y, gmcmcout, nmcmcout, cl_gmcmcout, cl_nmcmcout, prefix, group, subspl) {
  distrstats <- plot_diag_2distr(y, gmcmcout = cl_gmcmcout, nmcmcout = cl_nmcmcout,
                                 title = paste0(group, ' sample: ', subspl * 100, '%'))
  pdf(sprintf('%s_diag.pdf', prefix), width = 10, height = 3, onefile = FALSE)
  print(distrstats)
  dev.off()
  
  gamsum <- ggmcmcsummary(gmcmcout, cl_gmcmcout,
                          paste0(group, ' Gamma, sample: ', subspl * 100, '%'))
  pdf(sprintf('%s_gam_plots_gg.pdf', prefix), width = 7.8, height = 8.39, onefile = FALSE)
  print(gamsum)
  dev.off()
  
  norsum <- ggmcmcsummary(nmcmcout, cl_nmcmcout,
                          paste0(group, ' Normal, sample: ', subspl * 100, '%'))
  pdf(sprintf('%s_nor_plots_gg.pdf', prefix), width = 8.6, height = 8.6 / 2, onefile = FALSE)
  print(norsum)
  dev.off()
  
  gpac <- ggacplot(gmcmcout, paste0(group, ' Gamma full MCMC, sample: ', subspl * 100, '%'))
  cl_gpac <- ggacplot(cl_gmcmcout, paste0(group, ' Gamma burnin and thinning MCMC, sample: ', subspl * 100, '%'))
  pdf(sprintf('%s_gam_ac.pdf', prefix), width = 7, height = 7)
  print(gpac)
  print(cl_gpac)
  dev.off()

  npac <- ggacplot(nmcmcout, paste0(group, ' Normal full MCMC, sample: ', subspl * 100, '%'))
  cl_npac <- ggacplot(cl_nmcmcout, paste0(group, ' Normal burnin and thinning MCMC, sample: ', subspl * 100, '%'))
  pdf(sprintf('%s_nor_ac.pdf', prefix), width = 7, height = 7 / 2)
  print(npac)
  print(cl_npac)
  dev.off()
}
