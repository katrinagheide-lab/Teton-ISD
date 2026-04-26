#03 summaries of body size data

library(tidyverse)

# make sure path to your data is correct
dat <- read_rds("derived_data/TASR_dw.RDS")
nem <- read_rds("derived_data/TASR_nemouridae_dw.RDS")
dat <- bind_rows(dat, nem)
dat <- as_tibble(dat)
dat

# taxonomic counts
dat |>
  group_by(Label) |>
  summarize(n = n(), 
          min = min(dw),
          max = max(dw))

dat |>
  filter(dw > 10)

dat |>
  filter(Label == "Baetidae") |>
  distinct(file)

dat |>
  filter(Label == "Planaria") |>
  group_by(site, year) |>
  count()


# body size summaries
dat |>
  group_by(site, year) |>
  summarize(n = n(), 
            xmin = min(dw),
            xmedian = median(dw),
            xmean = mean(dw),
            xmax = max(dw),
            q10 = quantile(dw, probs = 0.1),
            q90 = quantile(dw, probs = 0.9))

ggplot(dat,
       aes(x = site, 
           y = dw, 
           color = as.factor(year))) +
  geom_boxplot()+
  scale_y_log10()

ggplot(dat,
       aes(x = site, 
           y = dw, 
           color = as.factor(year))) +
  geom_boxplot()+
  scale_y_log10() +
  facet_wrap(~Label)


ggplot(dat,
       aes(x = as.factor(year), 
           y = dw, 
           color = as.factor(year))) +
  geom_point(position = "jitter",
             alpha = 0.4,
             size = 1.5)+
  scale_y_log10() +
  facet_wrap(~site)

dat |>
  filter(dw < 0.015)
