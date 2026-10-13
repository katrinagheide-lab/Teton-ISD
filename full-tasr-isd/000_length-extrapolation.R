# 000 figuring out length extrapolations

# This script is to figure out how to extrapolate sub-sample measurements to full data counts
# Data has up to 20 individual measurements per taxon and sample
# There is often more than 20 individuals in a sample
# extrapolate 20 known measurements to full count data

library(tidyverse)

dat_2024 <- read_csv("annual_community_data/2024_inverts.csv")
dat_2024
dat_2024 |>
  distinct(Stream)
dat_2024 |>
  distinct(Subsample...8)
dat_2024 |>
  group_by(Stream) |>
  count()
dat_2024 |>
  distinct(Stage)
dat_2024 |>
  filter(Stage == "adult") |>
  distinct(Taxa)

cv <- dat_2024 |>
  filter(Stream == "Cloudveil Dome",
         !Stage %in% c("adult", "pupae"))
cv
names(cv)
cv <- cv |>
  select(Year:L20, Abundance:Density)
cv

cv |>
  mutate(non_na_count = rowSums(!is.na(pick(L1:L20)))) |>
  select(Taxa, Abundance, non_na_count) |>
  mutate(no_meas = Abundance - non_na_count,
         w = round(Abundance / non_na_count))

cv |>
  filter(Taxa == "Non-Tanypodinae") |>
  pivot_longer(L1:L20) |>
  group_by(Rep) |>
  reframe(quantile(value, probs = c(0.1, 0.5, 0.9)))

cv |>
  filter(Taxa == "Non-Tanypodinae") |>
  pivot_longer(L1:L20) |>
  ggplot(aes(x = value,
             fill = Rep)) +
  geom_histogram(position = "dodge",
                 binwidth = 1)

cv |>
  filter(Taxa == "Lednia") |>
  pivot_longer(L1:L20) |>
  filter(!is.na(value)) |>
  select(Taxa, Abundance, Rep, name, value) |>
  group_by(Rep, value) |>
  count()

cv |>
  filter(Taxa == "Lednia") |>
  pivot_longer(L1:L20) |>
  filter(!is.na(value)) |>
  select(Taxa, Abundance, Rep, name, value) |>
  group_by(Rep) |>
  mutate(n = n(),
         replace = Abundance>n,
         prop = Abundance / n) |>
  group_modify(~ slice_sample(.x,
                              prop = unique(.x$prop),
                              replace = unique(.x$replace))) |>
  group_by(Rep, value) |>
  count()


cv |>
  filter(Taxa == "Lednia") |>
  pivot_longer(L1:L20) |>
  filter(!is.na(value)) |>
  select(Taxa, Abundance, Rep, name, value) |>
  group_by(Rep) |>
  mutate(n = n(),
         replace = Abundance>n,
         prop = Abundance / n) |>
  group_modify(~ slice_sample(.x,
                              prop = unique(.x$prop),
                              replace = unique(.x$replace))) |>
  ggplot(aes(x = value,
             fill = Rep)) +
  geom_histogram(position = "dodge", 
                 binwidth = 0.5)

# all CV data
cv |>
  #filter(Taxa == "Lednia") |>
  pivot_longer(L1:L20) |>
  filter(!is.na(value)) |>
  select(Taxa, Abundance, Rep, name, value) |>
  group_by(Taxa, Rep) |>
  mutate(n = n(),
         replace = Abundance>n,
         prop = Abundance / n) |>
  group_modify(~ slice_sample(.x,
                              prop = unique(.x$prop),
                              replace = unique(.x$replace))) |>
  ggplot(aes(x = value,
             fill = Rep)) +
  geom_histogram(position = "dodge", 
                 binwidth = 0.5)

cv_sim <- cv |>
  pivot_longer(L1:L20) |>
  filter(!is.na(value)) |>
  select(Taxa, Abundance, Rep, name, value) |>
  group_by(Taxa, Rep) |>
  mutate(n = n(),
         replace = Abundance>n,
         prop = Abundance / n) |>
  group_modify(~ slice_sample(.x,
                              prop = unique(.x$prop),
                              replace = unique(.x$replace)))
cv_sim

# full data
dim(dat_2024)
# 300, 58
dat_2024 |>
  pivot_longer(L1:L20) |>
  filter(!is.na(value)) |>
  select(Stream, Fraction, Taxa, Abundance, Rep, name, value) |>
  group_by(Stream, Taxa, Rep) |>
  mutate(n = n(),
         replace = Abundance>n,
         prop = Abundance / n) |>
  View()
# middle teton rep 3 - duplicated rows?
dat_2024 |> group_by(Stream) |> distinct(Rep) |> arrange(Stream, Rep) |> View()

# how to figure out duplicates?
dat_2024 |>
  filter(!Stage %in% c("adult", "pupae")) |>
  group_by(Stream, Rep, Taxa) |>
  count() |>
  filter(n > 1)

  
# need to filter out the stream below
# otherwise, there are multiple "unique" combinations of stream, taxa, prop below
# function doesn't know how to deal with duplicates
dat_2024_sim <- dat_2024 |>
  filter(Stream != "Middle Teton",
         !Stage %in% c("adult", "pupae")) |>
  pivot_longer(L1:L20) |>
  filter(!is.na(value)) |>
  select(Stream, Fraction, Taxa, Abundance, Rep, name, value) |>
  group_by(Stream, Taxa, Rep) |>
  mutate(n = n(),
         replace = Abundance>n,
         prop = Abundance / n) |>
  group_modify(~ slice_sample(.x,
                              prop = unique(.x$prop),
                              replace = unique(.x$replace)))

dat_2024_sim|>
  ggplot(aes(x = value,
             fill = Rep)) +
  geom_histogram(position = "dodge", 
                 binwidth = 0.5) +
  facet_wrap(~Stream,
             scales = "free")

dat_2023 <- read_csv("annual_community_data/2023_inverts.csv")
dat_2023 |>
  filter(!Stage %in% c("adult", "pupae")) |>
  group_by(Stream, Rep, Fraction, Taxa) |>
  count() |>
  filter(n > 1)
# The Gusher rep = 1 Small --> duplicated rows?


# Cloudveil Dome Small rep 2 Lednia doeesn't have a total abundance
dat_2023_sim <- dat_2023 |>
  mutate(Abundance = case_when(
    Stream == "Cloudveil Dome" & Fraction == "Small" & Rep == 2 ~7, 
    .default = Abundance
  )) |>
  filter(Stream != "The Gusher",
         !Stage %in% c("adult", "pupae")) |>
  pivot_longer(L1:L20) |>
  filter(!is.na(value)) |>
  select(Stream, Fraction, Taxa, Abundance, Rep, name, value) |>
  group_by(Stream, Taxa, Rep, Fraction) |>
  mutate(n = n(),
         replace = Abundance>n,
         prop = Abundance / n) |>
  group_modify(~ slice_sample(.x,
                              prop = unique(.x$prop),
                              replace = unique(.x$replace)))

dat_2023_sim|>
  ggplot(aes(x = value,
             fill = as.character(Rep))) +
  geom_histogram(position = "dodge", 
                 binwidth = 0.5) +
  facet_wrap(~Stream,
             scales = "free")
