#02b_MLEbin

library(sizeSpectra)
# install.packages("poweRlaw")
# library(poweRlaw)
library(tidyverse)

dat <- read_rds("derived_data/TASR_dw.RDS")
nem <- read_rds("derived_data/TASR_nemouridae_dw.RDS")
dat <- bind_rows(dat, nem)
dat <- as_tibble(dat)
dat

dat |> 
  group_by(site, year) |>
  count()


mle_dat_list <- dat |>
  filter(dw>=0.0387) |>
  group_by(site, year) |>
  group_split()


# set up the for loop and data frame to store results
n = length(mle_dat_list)
n

# data frame results
bin_result <- data.frame(
  site = character(n),
  year = numeric(n),
  b = numeric(n), 
  CI_low = numeric(n),
  CI_high = numeric(n),
  xmin = numeric(n),
  xmax = numeric(n),
  n_x = numeric(n))

bin_result

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
  
  dw_bin <-  binData(x = body_mass,
                       binWidth = "2k")
  
  num.bins <- nrow(dw_bin$binVals)
  binBreaks <- c((dw_bin$binVals$binMin),
                 (dw_bin$binVals$binMax)[num.bins])
  binCounts <- dw_bin$binVals$binCount
  
  
  mle_bin <- calcLike(
    negLL.fn = negLL.PLB.binned,
    p = -2.5,
    w = binBreaks,
    d = binCounts,
    J = num.bins,
    vecDiff = 5,
    suppress.warnings = TRUE)
  
  mle_bin 
  
  
  
  # mle_estimate <- try(calcLike(
  #   negLL.fn = negLL.PLB,
  #   vecDiff = 5,
  #   x = body_mass, 
  #   xmin = xmin, 
  #   xmax = xmax, 
  #   n = n_x, 
  #   sumlogx = sum(log(body_mass)), 
  #   p = -1.5,
  #   suppress.warnings = TRUE))
  
  if(class(mle_bin) == "try-error"){
    bin_result$site[i] = site
    bin_result$year[i] = year
    bin_result$b[i] = NA
    bin_result$CI_low[i] = NA
    bin_result$CI_high[i] = NA
    bin_result$xmin[i] = NA
    bin_result$xmax[i] = NA
    bin_result$n_x[i] = NA
    next
  }
  
  if(exists("mle_bin")){
    # now save all the relevant information in the data frame
    # Make sure to put them in the right spot
    bin_result$site[i] = site
    bin_result$year[i] = year
    bin_result$b[i] = mle_bin$MLE
    bin_result$CI_low[i] = mle_bin$conf[1]
    bin_result$CI_high[i] = mle_bin$conf[2]
    bin_result$xmin[i] = xmin
    bin_result$xmax[i] = xmax
    bin_result$n_x[i] = n_x
    
    rm(mle_estimate)
  } 
}
bin_result



# plots -------------------------------------------------------------------

bin_result |>
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


site_mat <- tibble(site = rep(c("DELT", "GRIZ", "SCAS"), each = 2),
       year = rep(c(2018, 2024), 3),
       mat = c(1.75, 1.91,
               7.10, 10.71,
               2.15, 1.90))

bin_result <- bin_result |>
  left_join(site_mat)
  
bin_result |>
  ggplot(aes(x = mat, 
             y = b,
             ymin = CI_low, 
             ymax = CI_high,
             color = site)) +
  geom_pointrange(
    size = 1) +
  # scale_color_manual(values = c("#019AFF", "#FF1984")) +
  # scale_fill_manual(values = c("#019AFF", "#FF1984")) +
  theme_bw() +
  labs(y = "Estimated \U03BB") +
  #geom_line(position = position_dodge(width = 0.5)) +
  stat_smooth(inherit.aes = FALSE,
              aes(x = mat, y = b), 
              method = "lm") +
  NULL
  
summary(lm(b~mat, dat = bin_result))
summary(lme4::lmer(b~mat + (1|site) + (1|year), dat = bin_result))
