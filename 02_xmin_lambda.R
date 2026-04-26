# ruth powell 
# KAtrina Heide TASR

library(sizeSpectra)
# install.packages("poweRlaw")
library(poweRlaw)
library(tidyverse)

# custom function
# custom function 
MLE_tidy <- function(df, rsp_var){
  # define variables
  x <- df[[rsp_var]]
  x_n <- length(x)
  xmin = min(x)
  xmax = max(x)
  log.x = log(x)
  sum.log.x = sum(log.x)
  
  # initial starting point for parameter estimate
  PL.bMLE = 1/(log(min(x)) - sum.log.x/length(x)) - 1
  
  # non-linear minimization  
  PLB.minLL = nlm(negLL.PLB, 
                  p = PL.bMLE, x = x, n = length(x), 
                  xmin = xmin, xmax = xmax,
                  sumlogx = sum.log.x)
  
  # estimate for b
  PLB.bMLE = PLB.minLL$estimate
  # minimum estimate of b
  PLB.minNegLL = PLB.minLL$minimum
  
  ## 95% CI calculation
  bvec = seq(PLB.bMLE - 0.5, PLB.bMLE + 0.5, 1e-05)
  PLB.LLvals = vector(length = length(bvec))
  for (i in 1:length(bvec)) {
    PLB.LLvals[i] = negLL.PLB(bvec[i],
                              x = x,
                              n = length(x), 
                              xmin = xmin,
                              xmax = xmax,
                              sumlogx = sum.log.x)
  }
  critVal = PLB.minNegLL + qchisq(0.95, 1)/2
  bIn95 = bvec[PLB.LLvals < critVal]
  # confidence interval
  PLB.MLE.bConf = c(min(bIn95), max(bIn95))
  if (PLB.MLE.bConf[1] == min(bvec) | 
      PLB.MLE.bConf[2] == max(bvec)) {
    dev.new()
    plot(bvec, PLB.LLvals)
    abline(h = critVal, col = "red")
    stop("Need to make bvec larger - see R window")
  }
  # return b estimate and min/max 95% CI
  return(data.frame(b = PLB.bMLE,
                    minCI = min(bIn95),
                    maxCI = max(bIn95),
                    x_n = x_n))
}

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

# old mle estimate ####
# # actually estimate the isd exponent
# mle_lambda <- dat %>%
#   filter( dw >= xmin_clauset) |>
#   # filter out griz
#   filter(site != "GRIZ") |>
#   # mutate(keep = case_when(site == "GRIZ"& year == 2024 ~ 0 ,
#   #        .default = 1)) |>
#   # filter(keep == 1) |>
#   ## end of delet
#   group_by(site, year) %>% # site_number # NOT rep/sample
#   nest() %>%
#   mutate(lambda = map(data, # lambda = exponent/slope
#                       MLE_tidy,
#                       "dw")) %>%
#   unnest(cols = lambda) %>%
#   select(-data) %>%
#   ungroup()
# mle_lambda

# mle estimate with a list and loop ####
# based on tutorial at: https://jpomz.github.io/UFSCarISD/articles/tutorial.html#working-with-data-from-multiple-sites

mle_dat_list <- dat |>
  filter(dw >= 0.0387) |>
  #filter(dw >= xmin_clauset) |>
  group_by(site, year) %>% 
  group_split()
  
lapply(mle_dat_list, head)

# set up the for loop and data frame to store results
n = length(mle_dat_list)
n

# data frame results
isd_result <- data.frame(
  site = character(n),
  year = numeric(n),
  b = numeric(n), 
  CI_low = numeric(n),
  CI_high = numeric(n),
  xmin = numeric(n),
  xmax = numeric(n),
  n_x = numeric(n))

isd_result

# for loop
for(i in 1:n){
  # i = 4
  # read one piece of the list in at a time
  dat_in <- mle_dat_list[[i]] 
  # extract the relevant information from the dat_in object 
  body_mass <- dat_in$dw
  site <- unique(dat_in$site)
  year <- unique(dat_in$year)
  xmin = min(body_mass) 
  xmax = max(body_mass) 
  n_x = length(body_mass) 
  
  mle_estimate <- try(calcLike(
    negLL.fn = negLL.PLB,
    vecDiff = 5,
    x = body_mass, 
    xmin = xmin, 
    xmax = xmax, 
    n = n_x, 
    sumlogx = sum(log(body_mass)), 
    p = -1.5,
    suppress.warnings = TRUE))
  
  if(class(mle_estimate) == "try-error"){
    isd_result$site[i] = site
    isd_result$year[i] = year
    isd_result$b[i] = NA
    isd_result$CI_low[i] = NA
    isd_result$CI_high[i] = NA
    isd_result$xmin[i] = NA
    isd_result$xmax[i] = NA
    isd_result$n_x[i] = NA
    next
  }
  
  if(exists("mle_estimate")){
    # now save all the relevant information in the data frame
    # Make sure to put them in the right spot
    isd_result$site[i] = site
    isd_result$year[i] = year
    isd_result$b[i] = mle_estimate$MLE
    isd_result$CI_low[i] = mle_estimate$conf[1]
    isd_result$CI_high[i] = mle_estimate$conf[2]
    isd_result$xmin[i] = xmin
    isd_result$xmax[i] = xmax
    isd_result$n_x[i] = n_x
    
    rm(mle_estimate)
  } 
}
isd_result



# preliminary plots -------------------------------------------------------

isd_result |>
  #filter(site != "SCAS") |>
  ggplot(aes(x = year, 
             y = b,
             ymin = CI_low, 
             ymax = CI_high,
             color = site)) +
  geom_pointrange(
    size = 1,
    position = position_dodge(width = 0.5)) +
  # scale_color_manual(values = c("#019AFF", "#FF1984")) +
  # scale_fill_manual(values = c("#019AFF", "#FF1984")) +
  theme_bw() +
  labs(y = "Estimated \U03BB") +
  geom_line(position = position_dodge(width = 0.5)) +
  NULL



dat |>
  mutate(fill = dw >= xmin_clauset) |>
  #filter(dw < 7) |>
  ggplot(aes(x = (dw),
             fill = fill)) +
  geom_histogram(binwidth = 0.025) +
  facet_wrap(site~year,
             scales = "free_x") +
  #scale_y_log10() +
  scale_x_log10() +
  NULL

dat |>
  mutate(fill = dw >= 0.0387) |>
  #filter(dw < 7) |>
  ggplot(aes(x = (dw),
             fill = fill)) +
  geom_histogram(binwidth = 0.025) +
  facet_wrap(site~year,
             scales = "free_x") +
  #scale_y_log10() +
  scale_x_log10() +
  NULL


# without planaria?
dat |>
  mutate(fill = dw >= 0.0387) |>
  filter(Label != "Planaria",
         site == "GRIZ") |>
  ggplot(aes(x = (dw),
             fill = fill)) +
  geom_histogram(binwidth = 0.025) +
  facet_wrap(site~year,
             scales = "free_x") +
  #scale_y_log10() +
  scale_x_log10() +
  labs(title = "NO Planaria") +
  NULL

# ## old - update if lambda estimates improve
# x_01 <- seq(from = 0.196,
#             to = 76.7,
#             length.out = 806)
# line_01 = (1 - pPLB(x = x_01,
#                     b = -1.83, xmin = 0.196, xmax = 76.7)) * 806
# 
# x_02 <- seq(from = 0.196,
#             to = 6.86,
#             length.out = 806)
# line_02 = (1 - pPLB(x = x_02,
#                     b = -2.14, xmin = 0.196, xmax = 6.86)) * 806
# 
# line_01_df <- data.frame(x = x_01,
#                          y = line_01,
#                          line = "01")
# line_02_df <- data.frame(x = x_02,
#                          y = line_02,
#                          line = "03")
# line_dat <- bind_rows(line_01_df, line_02_df)
#   
# mle_dat %>% 
#   filter(dw >= 0.196) |>
#   group_by(site_number) %>%
#   sample_n(1000, replace = TRUE) |>
#   arrange(site_number, -dw) %>%
#   mutate(order = row_number()) %>%
#   ggplot(aes(x = dw,
#              y = order,
#              #size = body_mass,
#              color = site_number)) +
#   geom_point(shape = 1, alpha = 0.5) +
#   #theme_dark() +
#   theme_bw() +
#   scale_x_log10() +
#   scale_y_log10() +
#   labs(y = "Number of values \u2265 x",
#        x = "Individual body mass") +
#   geom_line(data = line_dat, 
#             aes(x = x, 
#                 y = y, 
#                 color = line),
#             inherit.aes = FALSE,
#             linewidth = 3, 
#             alpha = 0.5) +
#   ### add CI numbers here ###
#   annotate(geom = "text", x = 10, y = 0.8, label = "CI") +
#   ### add another annotate with other CI numbers here ###
#   NULL


# figure out log-log plots here ####
ggplot(dat,
       aes(x = dw, 
           fill = as.factor(year))) +
  geom_histogram(position = "identity",
                 alpha = 0.5) +
  scale_x_continuous(trans = "log2") +
  facet_wrap(~site)
bin_and_center <- function(data, var, breaks, ...){
  # data is a data frame
  # var is a string, and is the name of a column in data which you want to bin
  # breaks controls the number of bins as defined in hist() 
  # See ?hist for details
  
  # bin values using hist()
  binned_hist = hist(data[[var]], 
                     breaks = breaks, # need to predefine breaks
                     # e.g. Log2 breaks = 2^seq(min, max) 
                     # Log10 breaks = 10^seq(min, max)
                     include.lowest = TRUE, plot = FALSE)
  # calculate "left" and "right" edge of bins
  breaks_orig = binned_hist$breaks[1:(length(breaks)-1)]
  breaks_offset = binned_hist$breaks[2:length(breaks)]
  # total bin width = right edge - left edge
  break_width = breaks_offset - breaks_orig
  count = binned_hist$counts 
  dataout = data.frame(
    # normalize counts =count / width (White et al 1997)
    log_count_corrected = log10(count / break_width),
    # original midpoint of bin log10 transformed
    log_mids = log10(binned_hist$mids),
    log_mids_center = NA)
  # remove bins with 0 counts
  # -Inf comes from log10(count / break_width) above
  dataout = dataout[dataout$log_count_corrected !=-Inf,]
  # recenter data at x=0
  mid_row = ceiling(nrow(dataout)/2)
  # subtract value of mid row from all mids
  dataout$log_mids_center = 
    dataout[,"log_mids"] - dataout[mid_row,"log_mids"]
  dataout
}

breaks <- 2^seq(floor(log2(min(dat$dw))),
                ceiling(log2(max(dat$dw))))
breaks

log_lambda <- dat %>%
  filter( dw >= xmin_clauset) |>
  group_by(site, year) %>% # site_number # NOT rep/sample
  nest() %>%
  mutate(lambda = map(data, # lambda = exponent/slope
                      bin_and_center,
                      var = "dw",
                      breaks = breaks)) %>%
  unnest(cols = lambda) %>%
  select(-data) %>%
  ungroup() 

# log_lambda |>
#   ggplot(aes(x = log_mids, 
#              y = log_count_corrected,
#              color = as.factor(year),
#              shape = site,
#              fill = log_mids,
#              group = interaction(site, year))) +
#   geom_point(size = 5,
#              stroke = NA) +
#   scale_fill_viridis_c(option = "turbo") +
#   stat_smooth(aes(x = log_mids, 
#                   y = log_count_corrected,
#                   color = site),
#               inherit.aes = FALSE,
#               method = "lm",
#               se = FALSE) +
#   # scale_shape_manual(values = c(21, 24)) +
#   # scale_color_manual(values = c("red", "black")) +
#   theme_bw() + 
#   labs(y = "Log(Count)",
#        x = "Log(Body Size)") +
#   guides(fill = "none")


log_lambda |>
  ggplot(aes(x = log_mids, 
             y = log_count_corrected,
             color = as.factor(year))) +
  geom_point(size = 5,
             stroke = NA) +
  stat_smooth(method = "lm",
              se = FALSE) +
  # scale_shape_manual(values = c(21, 24)) +
  # scale_color_manual(values = c("red", "black")) +
  theme_bw() + 
  labs(y = "Log(Count)",
       x = "Log(Body Size)") +
  facet_wrap(~site)


log_lambda <- log_lambda |>
  group_by(site, year) |>
  mutate(id = cur_group_id())

log_lambda |>
  distinct(site, year, id)


summary(lm(log_count_corrected ~ 0+ log_mids_center * as.factor(id), data = log_lambda))
