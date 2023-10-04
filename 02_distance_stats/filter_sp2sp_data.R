get_spto <- function(x) {
  x[which(x != 'HUMAN')]
}

dat <- read.csv('../../01_distances/outputs_sp2sp/mammal_dist_sp2sp.csv')

datc <- dat[which(dat$MRCA_type == 'S' &
                    (dat$from_seq == 'HUMAN' | dat$to_seq == 'HUMAN') &
                    (dat$from == dat$tree | dat$to == dat$tree)), ]

datc[, 'sp_to'] <- apply(datc[, c(3, 5)], 1, get_spto)

datc <- datc[which(datc$ndist <= quantile(datc$ndist, 0.9)), ]

save(datc, file = '../data/sp2sp_dat.RData')
write.csv(datc, file = '../data/sp2sp_dat.csv', row.names = FALSE)
