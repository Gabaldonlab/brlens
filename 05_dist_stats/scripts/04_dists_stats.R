# Data analysis of phylome distances summary and relationships
# Moisès Bernabeu
# Barcelona, February, 2022

# Libraries loading
library(ggplot2)
theme_set(theme_bw())

# Data import
phy5 <- read.csv('../outputs/0005_stats_mrcasp.csv')
phy76 <- read.csv('../outputs/0076_stats_mrcasp.csv')

str(phy5)

# Data plot
ggplot(phy5, aes(Mean.dist, Mean.sp)) +
  geom_point()

ggplot(phy5, aes(Median.dist, Median.sp)) +
  geom_point()

ggplot(phy5, aes(Mean.dist, Mean.dupl)) +
  geom_point()

ggplot(phy5, aes(Median.dist, Median.dupl)) +
  geom_point()

ggplot(phy5, aes(Mean.dist_norm, Mean.sp)) +
  geom_point()

ggplot(phy5, aes(Median.dist_norm, Median.sp)) +
  geom_point()

ggplot(phy5, aes(Mean.dist_norm, Mean.dupl)) +
  geom_point()

ggplot(phy5, aes(Median.dist_norm, Median.dupl)) +
  geom_point()

ggplot(phy5, aes(Mean.dist, Mean.dupl / Mean.sp)) +
  geom_point()

ggplot(phy5, aes(Median.dist, Median.dupl / Median.sp)) +
  geom_point()

ggplot(phy76, aes(Mean.dist, Mean.sp)) +
  geom_point()

ggplot(phy76, aes(Median.dist, Median.sp)) +
  geom_point()

ggplot(phy76, aes(Mean.dist, Mean.dupl)) +
  geom_point()

ggplot(phy76, aes(Median.dist, Median.dupl)) +
  geom_point()

ggplot(phy76, aes(Mean.dist_norm, Mean.sp)) +
  geom_point()

ggplot(phy76, aes(Median.dist_norm, Median.sp)) +
  geom_point()

ggplot(phy76, aes(Mean.dist_norm, Mean.dupl)) +
  geom_point()

ggplot(phy76, aes(Median.dist_norm, Median.dupl)) +
  geom_point()

ggplot(phy76, aes(Mean.dist, Mean.dupl / Mean.sp)) +
  geom_point()

ggplot(phy76, aes(Median.dist, Median.dupl / Median.sp)) +
  geom_point()
