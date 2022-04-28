# Organising output data
# Moisès Bernabeu
# Barcelona, April 2022

# Definitions ----
get_other <- function(x, ref) {
  if (x['from_sp'] != ref) {
    y <- x['from_sp']
  } else {
    y <- x['to_sp']
  }
  return(y)
}

# Import and save image of data ----

# Read raw data
ryedat <- read.csv('../../02_get_distances/outputs/0005_dist.csv')
rhudat <- read.csv('../../02_get_distances/outputs/0076_dist.csv')

yespt <- read.csv('../../02_get_distances/outputs/0005_sptree_dist.csv')
huspt <- read.csv('../../02_get_distances/outputs/0076_sptree_dist.csv')

# Get species to
yedat <- ryedat
yedat$sp_to <- apply(ryedat, 1, get_other, ref = 'YEAST')
hudat <- rhudat
hudat$sp_to <- apply(rhudat, 1, get_other, ref = 'HUMAN')

yespt$sp_to <- apply(yespt, 1, get_other, ref = 'YEAST')
huspt$sp_to <- apply(huspt, 1, get_other, ref = 'HUMAN')

# Filter data
yedat <- yedat[which(yedat$mrca_type == 'S' & yedat$sp_to != 'YEAST' &
                       (yedat$to_sp == 'YEAST' | yedat$from_sp == 'YEAST') &
                       yedat$dist != 0), ]
hudat <- hudat[which(hudat$mrca_type == 'S' & hudat$sp_to != 'HUMAN' &
                       (hudat$to_sp == 'HUMAN' | hudat$from_sp == 'HUMAN') &
                       hudat$dist != 0), ]

yespt <- yespt[which(yespt$mrca_type == 'S' & yespt$sp_to != 'YEAST' &
                       (yespt$sp_to == 'YEAST' | yespt$from_sp == 'YEAST')), ]
huspt <- huspt[which(huspt$mrca_type == 'S' & huspt$sp_to != 'HUMAN' &
                       (huspt$sp_to == 'HUMAN' | huspt$from_sp == 'HUMAN')), ]


save(hudat, yedat, yespt, huspt, file = '../data/seed2sp_dist.Rdata')