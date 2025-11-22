devtools::install_github("limj0987/Tspecies")

library(Tspecies)

# The mutation rate is obtained from the paper: https://doi.org/10.1038/s41477-018-0323-6
miu <- 0.00000000302 # The unit is per site per year
# Conversions are shown below:
miu_g20 <- miu*20 # assume the generation time is 20 years
miu_g30 <- miu*30 # assume the generation time is 30 years

NACE <- readLines("NACE_ks_ind.ks")
NACE <- as.numeric(NACE)
cat('For L. tulipifera from North America (NA) and L. chinense from eastern China (CE): ')
Tspecies(NACE, miu_g20, 20)
Tspecies(NACE, miu_g30, 30)


NACW <- readLines("NACW_ks_ind.ks")
NACW <- as.numeric(NACW)
cat('For L. tulipifera from North America (NA) and L. chinense from western China (CW): ')
Tspecies(NACW, miu_g20, 20)
Tspecies(NACW, miu_g30, 30)

