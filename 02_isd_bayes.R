# 02_isd bayes


# install.packages("poweRlaw")
library(poweRlaw)
library(tidyverse)
library(brms)
library(isdbayes)
library(tidybayes)

# make sure path to your data is correct
dat <- read_rds("derived_data/TASR_dw.RDS")
nem <- read_rds("derived_data/TASR_nemouridae_dw.RDS")
dat <- bind_rows(dat, nem)
dat <- as_tibble(dat)
dat

##=================
# read in zapada data
# bind_rows etc. 


dat_list <- dat %>% 
  filter(!is.na(dw)) |>
  group_by(site, year) %>% 
  group_split()

# 3) create empty list to fill with xmins
xmin_list = list() 

# 4) get list of xmins for each sample
set.seed(202002)
for(i in 1:length(dat_list)){
  powerlaw = conpl$new(dat_list[[i]]$dw) # get power law estimate from poweRlaw package
  xmin_list[[i]] = tibble(xmin_clauset = estimate_xmin(powerlaw)$xmin, # extract the xmin from the poweRlaw package
                          site = unique(dat_list[[i]]$site),
                          year = unique(dat_list[[i]]$year))
}

(xmins_clauset = bind_rows(xmin_list))

dat <- left_join(dat, xmins_clauset) |>
  filter(!is.na(dw))
dat

dat_xmin <- dat |>
  #filter(dw >= xmin_clauset) |>
  group_by(site, year) |>
  mutate(x = dw,
         counts = 1, 
         xmin = min(dw), 
         xmax = max(dw)) |>
  select(site, year, x, counts, xmin, xmax)


fit = brm(x | vreal(counts, xmin, xmax) ~ site*year, 
           data = dat_xmin,
           stanvars = stanvars,
           family = paretocounts(),
           chains = 1, iter = 1000)

posts = fit$data |> 
  distinct(site, year, xmin, xmax) |> 
  mutate(counts = 1) |> 
  add_epred_draws(fit, re_formula = NA) 

posts |> 
  ggplot(aes(x = year,
             y = .epred,
             fill = site)) + 
  stat_halfeye(scale = 0.2) +
  theme_bw()
