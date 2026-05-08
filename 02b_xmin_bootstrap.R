# boot strap estimate for xmin
# written by JPZ
# based on a script from J Wesner


library(poweRlaw)
library(tidyverse)
# this method follows Clauset et al. 2009 by estimating the minimum
# size for which the data follow a power law using by minimizing the K-S statistic.
# i.e., xmin is model-based. 
# The procedure uses bootstrapping to get the geometric mean of xmins and then uses that geometric mean as the cutoff.

# load data
dat <-  read_rds("derived_data/TASR_all_dw.rds") |> 
  filter(dw > 0) |> 
  group_by(site, year) |> 
  # filter(correct_taxon != "Daphnia") |> 
  mutate(xmin = min(dw),
         xmax = max(dw)) 

dat
dat_list = dat |> group_by(site, year) |> group_split()

# Run bootstrap xmin procedure; it takes ~5 minutes. Completed version is loaded below.
xmin_list <- list()
k <- 1  # running index for results

for (i in seq_along(dat_list)) {
  for (rep in 1:100) {
    
    dat_i <- sample_n(dat_list[[i]],
                      replace = T,
                      size = nrow(dat_list[[i]])) # resample data
    powerlaw <- conpl$new(dat_i$dw) # fit power law from poweRlaw package
    
    xmin_list[[k]] <- tibble(
      xmin_clauset = estimate_xmin(powerlaw)$xmin,
      site = unique(dat_i$site),
      year = unique(dat_i$year),
      rep = rep,
      dataset = i
    ) # store the data
    
    k <- k + 1 # update the counter. Rinse and repeate
  }
}

xmin_list
# bind estimates and calculate quantiles
xmins_clauset = bind_rows(xmin_list) |>
  group_by(site, year) |> 
  rename(xmin_clauset_raw = xmin_clauset) |> 
  reframe(xmin_clauset = exp(mean(log(xmin_clauset_raw))),
          xmin_clauset_lower = quantile(xmin_clauset_raw, probs = 0.1),
          xmin_clauset_upper = quantile(xmin_clauset_raw, probs = 0.9))
xmins_boot <- xmins_clauset

ggplot(xmins_clauset,
       aes(x = year,
           color = site,
           y = xmin_clauset,
           ymin = xmin_clauset_lower,
           ymax = xmin_clauset_upper)) +
  geom_pointrange(position = position_jitter(width = 0.05)) +
  scale_y_log10() +
  geom_point(inherit.aes = FALSE,
             data = bind_rows(xmin_list),
             aes(x = year,
                 color = site,
                 y = xmin_clauset),
             size = 1, 
             alpha = 0.25,
             position = position_jitter(width = 0.05)) +
  theme_bw() +
  facet_wrap(year~site,
             scales = "free_x")

ggplot(xmins_clauset,
       aes(x = year,
           color = site,
           y = xmin_clauset,
           ymin = xmin_clauset_lower,
           ymax = xmin_clauset_upper)) +
  geom_point(shape = 18,
             size = 3) +
  scale_y_log10() +
  geom_point(inherit.aes = FALSE,
             data = bind_rows(xmin_list),
             aes(x = year,
                 color = site,
                 y = xmin_clauset),
             shape = 20,
             size = 1, 
             alpha = 0.25,
             position = position_jitter(width = 0.05)) +
  theme_bw() +
  facet_wrap(year~site,
             scales = "free_x")


