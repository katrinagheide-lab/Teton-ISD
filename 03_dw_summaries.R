#03 summaries of body size data

library(tidyverse)

# make sure path to your data is correct
dat <-  read_rds("derived_data/TASR_all_dw.rds")

# taxonomic counts
dat |>
  group_by(Label, year) |>
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

# predator biomass
distinct(dat, Label)
ffg_pred <- tibble(Label = c("Perlodidae",
                            "Rhyacophilidae",
                            "Simuliidae"),
                   pred = 1)


dat |>left_join(ffg_pred) |>
  mutate(pred = case_when(is.na(pred) ~ 0, .default = pred)) |>
  group_by(site, year, Label, pred) |>
  summarize(n = n(), 
            min = min(dw),
            mean = mean(dw),
            max = max(dw)) |> 
  arrange(-pred)

# total predator biomass
dat |> 
  left_join(ffg_pred) |>
  mutate(pred = case_when(is.na(pred) ~ 0, .default = pred)) |>
  group_by(site, year, pred) |>
  summarize(n = n(), 
            xmin = min(dw),
            xmedian = median(dw),
            xmean = mean(dw),
            xmax = max(dw),
            q10 = quantile(dw, probs = 0.1),
            q90 = quantile(dw, probs = 0.9)) |>
  filter(pred == 1) |>
  ggplot(aes(x = year, 
             y = xmean,
             ymin = q10, 
             ymax = q90,
             fill = site)) +
  geom_col(position = "dodge") +
  theme_bw() +
  labs(title = "Total Biomass of Predator Taxa",
       y = "Dry weight")


# proportion of predator biomass
# total predator biomass
dat |> 
  left_join(ffg_pred) |>
  mutate(pred = case_when(is.na(pred) ~ 0, .default = pred)) |>
  group_by(site, year, pred) |>
  summarize(biomass = sum(dw)) |>
  ungroup() |>
  group_by(site, year) |>
  mutate(tot_biomass = sum(biomass),
         prop_biomass = biomass / tot_biomass) |>
  #filter(pred == 1) |>
  ggplot(aes(x = year, 
             y = prop_biomass,
             fill = as.factor(pred))) +
  geom_col() +
  theme_bw() +
  labs(title = "Proportion of Biomass",
       y = "Dry weight") +
  facet_wrap(~site)
