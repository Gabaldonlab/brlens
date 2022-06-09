library(stringr)
library(tidyr)
library(latex2exp)

hutables <- list.files('/data/moises/cluster/documents/brlens/09_gamma_inference/sp2sp_outputs/', pattern = 'hudat.*.csv', full.names = TRUE)

getsp <- function(filename) {
  str_split(str_split(filename, '/', simplify = TRUE)[, 10], '_', simplify = TRUE)[, 3]
}

hudat <- data.frame()
for (datfile in hutables) {
  hudat <- rbind(hudat, data.frame(read.csv(datfile), spto = getsp(datfile)))
}

yetables <- list.files('/data/moises/cluster/documents/brlens/09_gamma_inference/sp2sp_outputs/', pattern = 'yedat.*.csv', full.names = TRUE)

yedat <- data.frame()
for (datfile in yetables) {
  yedat <- rbind(yedat, data.frame(read.csv(datfile), spto = getsp(datfile)))
}

yedat <- data.frame(yedat, dataset = 'Yeast')
hudat <- data.frame(hudat, dataset = 'Human')
dat <- rbind(yedat, hudat)
names(dat)[1] <- 'stat'

# write.csv(dat, file = '../outputs/posterior_data.csv')

dat <- read.csv('../outputs/posterior_data.csv')

# pdf('../outputs/mode_vs_pars.pdf', width = 7.25, height = 8.5)
par(mfrow = c(3, 2))
plot(dat[which(dat$stat == 'mo' & dat$dataset == 'Yeast'), 'Mean'],
     dat[which(dat$stat == 'a' & dat$dataset == 'Yeast'), 'Mean'],
     xlab = 'Posterior distribution for the mode mean',
     ylab = TeX('Posterior distribution for the $\\alpha$ mean'),
     main = 'Yeast phylome', pch = 20)
plot(dat[which(dat$stat == 'mo' & dat$dataset == 'Human'), 'Mean'],
     dat[which(dat$stat == 'a' & dat$dataset == 'Human'), 'Mean'],
     xlab = 'Posterior distribution for the mode mean',
     ylab = TeX('Posterior distribution for the $\\alpha$ mean'),
     main = 'Yeast phylome', pch = 20)

plot(dat[which(dat$stat == 'mo' & dat$dataset == 'Yeast'), 'Mean'],
     dat[which(dat$stat == 'b' & dat$dataset == 'Yeast'), 'Mean'],
     xlab = 'Posterior distribution for the mode mean',
     ylab = TeX('Posterior distribution for the $\\beta$ mean'),
     main = 'Human phylome', pch = 20)
plot(dat[which(dat$stat == 'mo' & dat$dataset == 'Human'), 'Mean'],
     dat[which(dat$stat == 'b' & dat$dataset == 'Human'), 'Mean'],
     xlab = 'Posterior distribution for the mode mean',
     ylab = TeX('Posterior distribution for the $\\beta$ mean'),
     main = 'Human phylome', pch = 20)

plot(dat[which(dat$stat == 'mo' & dat$dataset == 'Yeast'), 'Mean'],
     dat[which(dat$stat == 'v' & dat$dataset == 'Yeast'), 'Mean'],
     xlab = 'Posterior distribution for the mode mean',
     ylab = 'Posterior distribution for the Var. mean',
     main = 'Yeast phylome', pch = 20)
plot(dat[which(dat$stat == 'mo' & dat$dataset == 'Human'), 'Mean'],
     dat[which(dat$stat == 'v' & dat$dataset == 'Human'), 'Mean'],
     ylim = c(0, 15),
     xlab = 'Posterior distribution for the mode mean',
     ylab = 'Posterior distribution for the Var. mean',
     main = 'Human phylome', pch = 20)
# dev.off()

sdat <- rbind('Effective size' = summary(dat$Effective.size[-which(is.na(dat$R))]),
              'R' = summary(dat$R[-which(is.na(dat$R))]))
