#04_community_data
library(tidyverse)

# tasr 2024
tasr_2024 <- read_csv("derived_data/2024GRTEStonefly.csv")
tasr_2024
distinct(tasr_2024, Site)


# tasr 2015-2018
tasr_2015_2018 <- read_csv("derived_data/TASR_inverts_2015-2018.csv") 
distinct(tasr_2015_2018, Stream, Year) |>
  filter(Year == 2018)

rp18 <- tasr_2015_2018|>
  filter(Stream %in% c("Delta Inlet Stream",
                       "Grizzly Inlet Stream",
                       "South Cascade"),
         Year == 2018,
         !str_detect(Taxa, "pup"),
         Taxa != "Staphylinidae") |>
  select(Stream, Year, Fraction, Rep, Taxa, Density, Abundance)

rp18
distinct(rp18, Taxa)

# filter out sites for Ruth Powell analysis
rp24 <- tasr_2024 |>
  filter(Site %in% c("Delta Lake Inlet", "Grizzly", "South Cascade RG")) |>
  select(Year, Site, Rep, Subsample...8, Taxa, Abundance, Density_ind_m2, Biomass_mg_m2)

rp24
distinct(rp24, Taxa)
# Predator FFG pulled from Jorgeson et al. 2024
# supplemented with Merrit
ffg_pred <- tibble(Taxa = c("2Megarcys",
                            "Prosimulium",
                            "Megarcys",
                            "Sweltsa", 
                            "Rhyacophila",
                            "Rhyacophila Alberta Group",
                            "Rhyacophila blarina", 
                            "Rhyacophila Verrula Group",
                            "Simuliidae"),
                   Pred = 1)

# Abundance
# plot proportion of taxa which are classified as predators
# Abundance
rp24 |>
  left_join(ffg_pred) |>
  mutate(Pred = case_when(is.na(Pred) ~ 0, .default = Pred)) |>
  group_by(Year, Site, Taxa, Pred) |>
  summarize(mean_n = mean(Abundance, na.rm = TRUE),
            mean_dens = mean(Density_ind_m2, na.rm = TRUE),
            sd_dens = sd(Density_ind_m2, na.rm = TRUE)) |>
  ungroup()|>
  group_by(Year, Site)|>
  mutate(tot_n = sum(mean_dens),
         taxa_prop = mean_dens / tot_n) |>
  ggplot(aes(x = Site, 
             y = taxa_prop,
             fill = as.factor(Pred))) +
  geom_col()


  
# predators as a proportion of total abundance
rp24_pred <- rp24 |>
  left_join(ffg_pred) |>
  mutate(Pred = case_when(is.na(Pred) ~ 0, .default = Pred)) |>
  group_by(Year, Site, Taxa, Pred) |>
  summarize(mean_n = mean(Abundance, na.rm = TRUE),
            mean_dens = mean(Density_ind_m2, na.rm = TRUE),
            sd_dens = sd(Density_ind_m2, na.rm = TRUE)) |>
  ungroup()|>
  group_by(Year, Site)|>
  mutate(tot_n = sum(mean_dens),
         taxa_prop = mean_dens / tot_n) |>
  group_by(Year, Site, Pred) |>
  summarise(prop_pred = sum(taxa_prop)) |>
  filter(Pred == 1)
rp24_pred  
# Grizz in 2024 ~ 38% predators by abundance
# Delta and SCAS ~ 0.5-0.7%



# 2018 --------------------------------------------------------------------

rp18 |>
  left_join(ffg_pred) |>
  mutate(Pred = case_when(is.na(Pred) ~ 0, .default = Pred)) |>
  group_by(Year, Stream, Taxa, Pred) |>
  summarize(mean_n = mean(Abundance, na.rm = TRUE),
            mean_dens = mean(Density, na.rm = TRUE),
            sd_dens = sd(Density, na.rm = TRUE)) |>
  ungroup()|>
  group_by(Year, Stream)|>
  mutate(tot_n = sum(mean_dens),
         taxa_prop = mean_dens / tot_n) |>
  ggplot(aes(x = Stream, 
             y = taxa_prop,
             fill = as.factor(Pred))) +
  geom_col()

rp18 |>
  left_join(ffg_pred) |>
  mutate(Pred = case_when(is.na(Pred) ~ 0, .default = Pred)) |>
  group_by(Year, Stream, Taxa, Pred) |>
  summarize(mean_n = mean(Abundance, na.rm = TRUE),
            mean_dens = mean(Density, na.rm = TRUE),
            sd_dens = sd(Density, na.rm = TRUE)) |>
  ungroup()|>
  group_by(Year, Stream)|>
  mutate(tot_n = sum(mean_dens),
         taxa_prop = mean_dens / tot_n) |>
  group_by(Year, Stream, Pred) |>
  summarise(sum(taxa_prop))

# Grizzly ~17% Predators
# Delta ~ 1 %
# SCAS = 0%



# other years -------------------------------------------------------------

rp_all <- tasr_2015_2018|>
  filter(Stream %in% c("Delta Inlet Stream",
                       "Grizzly Inlet Stream",
                       "South Cascade"),
         #Year == 2018,
         !str_detect(Taxa, "pup"),
         Taxa != "Staphylinidae") |>
  select(Stream, Year, Fraction, Rep, Taxa, Density, Abundance)

distinct(rp_all, Taxa)

ffg_pred <- tibble(Taxa = c("2Megarcys",
                            "Megarcys2",
                            "1Megarcys",
                            "Alloperla", #Chloroperlidae family generally predators
                            "Chloroperlidae",
                            "Perlodidae",
                            "Clinocera", # Empididae generally predators
                            "Suwallia",
                            "Drunella",
                            "Empididae",
                            "Megarcys",
                            "Sweltsa", 
                            "Rhyacophila",
                            "Rhyacophila Alberta Group",
                            "Rhyacophila blarina", 
                            "Rhyacophila Verrula Group",
                            "Prosimulium",
                            "Simulium",
                            "Simuliidae"),
                   Pred = 1)

rp_all |>
  left_join(ffg_pred) |>
  mutate(Pred = case_when(is.na(Pred) ~ 0, .default = Pred)) |>
  group_by(Year, Stream, Taxa, Pred) |>
  summarize(mean_n = mean(Abundance, na.rm = TRUE),
            mean_dens = mean(Density, na.rm = TRUE),
            sd_dens = sd(Density, na.rm = TRUE)) |>
  ungroup()|>
  group_by(Year, Stream)|>
  mutate(tot_n = sum(mean_dens),
         taxa_prop = mean_dens / tot_n) |>
  ggplot(aes(x = Year, 
             y = taxa_prop,
             fill = as.factor(Pred))) +
  geom_col() +
  facet_wrap(~Stream)

rp_all |>
  left_join(ffg_pred) |>
  mutate(Pred = case_when(is.na(Pred) ~ 0, .default = Pred)) |>
  group_by(Year, Stream, Taxa, Pred) |>
  summarize(mean_n = mean(Abundance, na.rm = TRUE),
            mean_dens = mean(Density, na.rm = TRUE),
            sd_dens = sd(Density, na.rm = TRUE)) |>
  ungroup()|>
  group_by(Year, Stream)|>
  mutate(tot_n = sum(mean_dens),
         taxa_prop = mean_dens / tot_n) |>
  group_by(Year, Stream, Pred) |>
  summarise(prop_pred = sum(taxa_prop)) |>
  filter(Pred == 1) |>
  mutate(Site = case_when(
    Stream == "South Cascade" ~ "South Cascade RG",
    Stream == "Delta Inlet Stream" ~ "Delta Lake Inlet",
    Stream == "Grizzly Inlet Stream" ~ "Grizzly")) |>
  bind_rows(rp24_pred) |>
  ggplot(aes(x = Year, 
             y = prop_pred,
             color = Site)) +
  geom_line() +
  geom_point()
