# 003 - length weight conversions

library(tidyverse)
# read in length data sets
dat1 <- readRDS("full-tasr-isd/cleaned_length_data/2015-2018.rds")
dat1$Year <- as.numeric(dat1$Year)
dat2 <- readRDS("full-tasr-isd/cleaned_length_data/2019-2024.rds")
# bind two data sets together
dat <- bind_rows(dat1, dat2)

# read in LW coeffs
lw_coef <- read.csv("LW_coeffs.csv")

lw_cols <- c("taxon",
             "a",
             "b",
             "L_units",
             "dw_units",
             "formula_type",
             "formula")

lw_coef <- lw_coef[, names(lw_coef) %in% lw_cols]
names(lw_coef)

# compare equations with taxa list
taxa_dat <- dat |> distinct(Taxa) |> pull() |> sort()
taxa_eq <- lw_coef |> distinct(taxon) |> pull() |> sort()
setdiff(taxa_dat, taxa_eq)

dat <- dat |>
  mutate(taxon = case_when(
    Taxa == "Alloperla" ~"Chloroperlidae",
    Taxa == "Agathon" ~"Blephariceridae",
    Taxa == "Allomyia" ~"Trichoptera",
    Taxa == "Bezzia" ~"Ceratopogonidae",
    Taxa == "Celina" ~"Coleoptera",
    Taxa == "Cinygmula" ~"Heptageniidae",
    Taxa == "Cleptelmis" ~"Elmidae",
    Taxa == "Clinocera" ~"Empididae",
    Taxa == "Culicoides" ~"Ceratopogonidae",
    Taxa == "Dicranota" ~"Tipulidae",
    Taxa == "Dixella" ~"Diptera",
    Taxa == "Drunella " ~"Ephemerellidae",
    Taxa == "Ecclisomyia" ~"Limnephilidae",
    Taxa == "Glossosoma" ~"Glossosomatidae",
    Taxa == "Gonomyodes" ~"Tipulidae",
    Taxa == "Helodon" ~"Simuliidae",
    Taxa == "Homophylax" ~"Limnephilidae",
    Taxa == "Skwala" ~"Perlodidae",
    Taxa == "Lednia" ~"Nemouridae",
    Taxa == "Limonia" ~"Tipulidae",
    Taxa == "Limoniidae" ~"Tipulidae",
    Taxa == "Rhyacophila" ~"Rhyacophilidae",
    Taxa == "Megarcys" ~"Perlodidae",
    Taxa == "Molophilus" ~"Tipulidae",
    Taxa == "Rhithrogena" ~"Heptageniidae",
    Taxa == "Oreogeton" ~"Empididae",
    Taxa == "Ormosia" ~"Tipulidae",
    Taxa == "Paradelphomyia" ~"Tipulidae",
    Taxa == "Parapsyche" ~"Hydropsychidae",
    Taxa == "Psychoglypha" ~"Limnephilidae",
    Taxa == "Pomoleuctra" ~"Leuctridae",
    Taxa == "Probezzia" ~"Ceratopogonidae",
    Taxa == "Prosimuliium" ~"Simuliidae",
    Taxa == "Stenelmis" ~"Elmidae",
    Taxa == "Suwallia" ~"Chloroperlidae",
    Taxa == "Sweltsa" ~"Chloroperlidae",
    Taxa == "Tanypodinae" ~"Chironomidae",
    Taxa == "Twinnia" ~"Simuliidae",
    Taxa == "Visoka" ~"Nemouridae",
    Taxa == "Zapada" ~"Nemouridae",
    Taxa == "Yoraperla" ~"Peltoperlidae",
    .default = Taxa
  ))
taxa_dat_corrected <- dat |> distinct(taxon) |> pull() |> sort()
no_eq <- setdiff(taxa_dat_corrected, taxa_eq)

dat |>
  filter(taxon %in% no_eq)
# ~206 rows without equations
206 / nrow(dat) *100
# ~ 0.18% without equations

dat_eq <- dat |>
  inner_join(lw_coef)
# estimate dw based on formula type ####
dat_mass <- dat_eq |>
  filter(l_mm > 2.5) |>
  mutate(dw = case_when(
    formula_type == 1 ~ a * l_mm^b,
    formula_type == 2 ~ exp(a + b * log(l_mm)))) |>
  select(Year, Stream, taxon, dw)
dat_mass
