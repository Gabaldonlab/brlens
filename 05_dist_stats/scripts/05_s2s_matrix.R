dat <- read.csv('../data/0005_dists_noh.csv')

sp2sp.dat <- dat[which(dat$mrca_type == 'S'), ]

sps <- dat$from_sp[!duplicated(dat$from_sp)]
sp.com <- t(combn(sps, 2))

sp2sp.mat <- array(dim = c(length(sps), length(sps), 4),
                   dimnames = list(sps, sps, c('dist.mea', 'dist.med',
                                               'ndist.mea', 'ndist.med')))
for (i in 1:dim(sp.com)[1]) {
  print(i)
  tr.dat <- dat[which(dat$from_sp == sp.com[i, 1] & dat$to_sp == sp.com[i, 2] &
                        dat$mrca_type == 'S'), ]
  sp2sp.mat[sp.com[i, 1], sp.com[i, 2], 'dist.mea'] <- mean(tr.dat$dist)
  sp2sp.mat[sp.com[i, 1], sp.com[i, 2], 'dist.med'] <- median(tr.dat$dist)
  sp2sp.mat[sp.com[i, 1], sp.com[i, 2], 'ndist.mea'] <- mean(tr.dat$dist_norm)
  sp2sp.mat[sp.com[i, 1], sp.com[i, 2], 'ndist.med'] <- median(tr.dat$dist_norm)
  sp2sp.mat[sp.com[i, 2], sp.com[i, 1], 'dist.mea'] <- mean(tr.dat$dist)
  sp2sp.mat[sp.com[i, 2], sp.com[i, 1], 'dist.med'] <- median(tr.dat$dist)
  sp2sp.mat[sp.com[i, 2], sp.com[i, 1], 'ndist.mea'] <- mean(tr.dat$dist_norm)
  sp2sp.mat[sp.com[i, 2], sp.com[i, 1], 'ndist.med'] <- median(tr.dat$dist_norm)
}

heatmap(sp2sp.mat[, , 'dist.mea'], main = 'dist.mea')
heatmap(sp2sp.mat[, , 'dist.med'], main = 'dist.med')
heatmap(sp2sp.mat[, , 'ndist.mea'], main = 'ndist.mea')
heatmap(sp2sp.mat[, , 'ndist.med'], main = 'ndist.med')
