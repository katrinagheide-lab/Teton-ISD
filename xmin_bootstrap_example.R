library(poweRlaw)
library(tidyverse)
# this method follows Clauset et al. 2009 by estimating the minimum
# size for which the data follow a power law using by minimizing the K-S statistic.
# i.e., xmin is model-based. 
# The procedure uses bootstrapping to get the geometric mean of xmins and then uses that geometric mean as the cutoff.

# load data
dat = read_csv("data/dw_fixed.csv") %>% 
  filter(dw_mg > 0) %>% 
  group_by(tank) %>% 
  # filter(correct_taxon != "Daphnia") %>% 
  mutate(xmin = min(dw_mg),
         xmax = max(dw_mg)) 

dat_list = dat %>% group_by(tank) %>% group_split()

# Run bootstrap xmin procedure; it takes ~5 minutes. Completed version is loaded below.
xmin_list <- list()
k <- 1  # running index for results

for (i in seq_along(dat_list)) {
  for (rep in 1:100) {

    dat_i <- sample_n(dat_list[[i]], replace = T, size = nrow(dat_list[[i]])) # resample data
    powerlaw <- conpl$new(dat_i$dw_mg) # fit power law from poweRlaw package

    xmin_list[[k]] <- tibble(
      xmin_clauset = estimate_xmin(powerlaw)$xmin,
      tank = unique(dat_i$tank),
      rep = rep,
      dataset = i
    ) # store the data

    k <- k + 1 # update the counter. Rinse and repeate
  }
}
saveRDS(xmin_list, file = "models/xmin_list.rds")

xmins_clauset = bind_rows(readRDS("models/xmin_list.rds")) %>% group_by(tank) %>% 
  rename(xmin_clauset_raw = xmin_clauset) %>% 
  reframe(xmin_clauset = exp(mean(log(xmin_clauset_raw))),
          xmin_clauset_lower = quantile(xmin_clauset_raw, probs = 0.1),
          xmin_clauset_upper = quantile(xmin_clauset_raw, probs = 0.9))

saveRDS(xmins_clauset, file = "models/xmin_boots.rds")
xmins_clauset = readRDS(file = "models/xmin_boots.rds")


# cull data to remove body sizes below xmin_clauset and recalculate xmin and xmax
# this creates the data used in model fitting in 04_fit_models.R
dat_clauset_xmins_boot = dat %>% left_join(xmins_clauset) %>% 
  group_by(tank) %>% 
  filter(dw_mg >= xmin_clauset) %>%
  mutate(xmin = min(dw_mg),
         xmax = max(dw_mg),
         source = "bootstrapped_xmins")

saveRDS(dat_clauset_xmins_boot, file = "data/dat_clauset_xmins.rds")

dat_clauset_xmins_boot_lower = dat %>% left_join(xmins_clauset) %>% 
  group_by(tank) %>% 
  filter(dw_mg >= xmin_clauset_lower) %>%
  mutate(xmin = min(dw_mg),
         xmax = max(dw_mg),
         source = "bootstrapped_xmins")

saveRDS(dat_clauset_xmins_boot_lower, file = "data/dat_clauset_xmins_lower.rds")

dat_clauset_xmins_boot_upper = dat %>% left_join(xmins_clauset) %>% 
  group_by(tank) %>% 
  filter(dw_mg >= xmin_clauset_upper) %>%
  mutate(xmin = min(dw_mg),
         xmax = max(dw_mg),
         source = "bootstrapped_xmins")

saveRDS(dat_clauset_xmins_boot_upper, file = "data/dat_clauset_xmins_upper.rds")

# Remove Daphnia ----------------------------------------------------------

library(poweRlaw)
library(tidyverse)
# this method follows Clauset et al. 2009 by estimating the minimum
# size for which the data follow a power law using by minimizing the K-S statistic.
# i.e., xmin is model-based. 
# The procedure uses bootstrapping to get the geometric mean of xmins and then uses that geometric mean as the cutoff.

# load data
dat_nodaphnia = read_csv("data/dw_fixed.csv") %>% 
  filter(dw_mg > 0) %>% 
  group_by(tank) %>% 
  filter(correct_taxon != "Daphnia") %>%
  mutate(xmin = min(dw_mg),
         xmax = max(dw_mg)) 

dat_list_nodaphnia = dat_nodaphnia %>% group_by(tank) %>% group_split()
# Run bootstrap xmin procedure; silenced b/c it takes ~5 minutes. Completed version is loaded below.
xmin_list_nodaphnia <- list()
k <- 1  # running index for results

for (i in seq_along(dat_list_nodaphnia)) {
  for (rep in 1:100) {
    
    dat_i <- sample_n(dat_list_nodaphnia[[i]], replace = T, size = nrow(dat_list_nodaphnia[[i]]))
    powerlaw <- conpl$new(dat_i$dw_mg)
    
    xmin_list_nodaphnia[[k]] <- tibble(
      xmin_clauset = estimate_xmin(powerlaw)$xmin,
      tank = unique(dat_i$tank),
      rep = rep,
      dataset = i
    )
    
    k <- k + 1
  }
}
saveRDS(xmin_list_nodaphnia, file = "models/xmin_list_nodaphnia.rds")

xmins_clauset_nodaphnia = bind_rows(readRDS("models/xmin_list_nodaphnia.rds")) %>% group_by(tank) %>% 
  rename(xmin_clauset_raw = xmin_clauset) %>% 
  reframe(xmin_clauset = exp(mean(log(xmin_clauset_raw))),
          xmin_clauset_lower = quantile(xmin_clauset_raw, probs = 0.1),
          xmin_clauset_upper = quantile(xmin_clauset_raw, probs = 0.9))

saveRDS(xmins_clauset_nodaphnia, file = "models/xmin_boots_nodaphnia.rds")
xmins_clauset_nodaphnia = readRDS(file = "models/xmin_boots_nodaphnia.rds")



# data with emergers --------------------
# emergence data collected 27-28 June at the end of the experiment. Same time as benthics
tictoc::tic()
dw_raw_adults = read_csv("data/dw_raw_adults.csv") %>% select(tank, temp_treat, nutrient_treat, dw_mg, m2) %>% 
  mutate(stage = "adult")
dat = read_csv("data/dw_fixed.csv") %>% filter(dw_mg > 0) %>% mutate(m2 = 0.09) %>% select(tank, temp_treat, nutrient_treat, dw_mg, m2, correct_taxon) %>% 
  mutate(stage = "larvae") %>% 
  filter(correct_taxon != "Daphnia")

dat_emerge_ind = bind_rows(dw_raw_adults, dat) %>% 
  group_by(tank, temp_treat, nutrient_treat, dw_mg, m2, stage) %>% 
  tally() %>% 
  mutate(no_m2 = n/m2)

dat_list_emerge = dat_emerge_ind %>% group_by(tank) %>% group_split()

xmin_list_emerge = list()

k <- 1  # running index for results

for (i in seq_along(dat_list_emerge)) {
  for (rep in 1:100) {

    dat_i <- sample_n(dat_list_emerge[[i]], replace = T, size = nrow(dat_list_emerge[[i]]))
    powerlaw <- conpl$new(dat_i$dw_mg)

    xmin_list_emerge[[k]] <- tibble(
      xmin_clauset = estimate_xmin(powerlaw)$xmin,
      tank = unique(dat_i$tank),
      rep = rep,
      dataset = i
    )

    k <- k + 1
  }
}

xmins_clauset_emerge = bind_rows(xmin_list_emerge) %>% group_by(tank) %>% reframe(xmin_clauset = exp(mean(log(xmin_clauset))))
saveRDS(xmins_clauset_emerge, file = "models/xmin_boots_emerge.rds")

dat_clauset_xmins_emerge = dat_emerge_ind %>% left_join(xmins_clauset_emerge) %>% 
  group_by(tank) %>% 
  filter(dw_mg >= xmin_clauset) %>%
  mutate(xmin = min(dw_mg),
         xmax = max(dw_mg))

saveRDS(dat_clauset_xmins_emerge, file = "data/dat_clauset_xmins_emerge.rds")
tictoc::toc()

