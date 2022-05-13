dat <- read.csv('../outputs/rand_phylome_dist.csv')

spdat <- dat[which(dat$from == 'A' | 
                     dat$to == 'A' & dat$from != dat$to), ]

mrca_ref <- read.csv('../outputs/mrca.csv', row.names = 1)
rbls_ref <- read.csv('../outputs/rbls.csv', row.names = 1)
root_ref <- read.csv('../outputs/root.csv', row.names = 1)
st_ref <- read.csv('../outputs/st.csv', row.names = 1)

spdat$mrca_dist <- spdat$dist / mrca_ref[spdat$id, 'med'] / median(mrca_ref[, 'med'])
spdat$root_dist <- spdat$dist / root_ref[spdat$id, 'med'] / median(root_ref[, 'med'])
spdat$st_dist <- spdat$dist / st_ref[spdat$id, 'med'] / median(st_ref[, 'med'])
spdat$rbls_sum <- spdat$dist / rbls_ref[spdat$id, 'sum_brl'] / median(rbls_ref[, 'sum_brl'])
spdat$rbls_med <- spdat$dist / rbls_ref[spdat$id, 'med_brl'] / median(rbls_ref[, 'med_brl'])

spdatcl <- spdat[, -c(1:2)]
str(spdat)

spdatcl <- gather(spdatcl, key = 'type', value = 'dist', -to)

ggplot(spdatcl, aes(x = dist, col = to)) +
  geom_density() +
  facet_wrap(~type, scales = 'free')

a <- ggplot(spdat, aes(x = dist, col = to)) +
  geom_density() +
  xlab('Raw distance') +
  ylab('Density') +
  labs(colour = 'Species to')

b <- ggplot(spdat, aes(x = st_dist, col = to)) +
  geom_density() +
  xlab('Subtree n. distance') +
  ylab('Density') +
  labs(colour = 'Species to') +
  xlim(0, quantile(spdat$st_dist, 0.9))

c <- ggplot(spdat, aes(x = mrca_dist, col = to)) +
  geom_density() +
  xlab('MRCA n. distance') +
  ylab('Density') +
  labs(colour = 'Species to') +
  xlim(0, quantile(spdat$st_dist, 0.95))

d <- ggplot(spdat, aes(x = root_dist, col = to)) +
  geom_density() +
  xlab('Root to tip n. distance') +
  ylab('Density') +
  labs(colour = 'Species to') +
  xlim(0, quantile(spdat$st_dist, 0.9))

e <- ggplot(spdat, aes(x = rbls_med, col = to)) +
  geom_density() +
  xlab('Br. len. med n. distance') +
  ylab('Density') +
  labs(colour = 'Species to') +
  xlim(0, quantile(spdat$st_dist, 0.999))

f <- ggplot(spdat, aes(x = rbls_sum, col = to)) +
  geom_density() +
  xlab('Br. len. sum n. distance') +
  ylab('Density') +
  labs(colour = 'Species to') +
  xlim(0, quantile(spdat$st_dist, 0.2))

# pdf('../outputs/phylome_norms.pdf', width = 9, height = 4.8)
ggarrange(a, b, c, d, e, f, common.legend = TRUE,
          legend = 'bottom', labels = 'auto')
# dev.off()
