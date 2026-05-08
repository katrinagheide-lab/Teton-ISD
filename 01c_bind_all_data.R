# 01c_stitch all data
library(tidyverse)

# read in most of the data
dat <- read_rds("derived_data/TASR_dw.RDS")
# read in the simulated nemouridae data
nem <- read_rds("derived_data/TASR_nemouridae_dw.RDS")
# bind them together
dat <- bind_rows(dat, nem)
dat <- as_tibble(dat)
dat

# save the "full" data set
saveRDS(dat, 
        "derived_data/TASR_all_dw.rds")
