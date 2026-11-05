# 2019-2024 data
# this script is looking to clean and prepare data for isdbayes analysis
# standardize data format
# look for duplicates and errors
# clean data
# pivot longer
# extrapolate lengths
# prepare data for LW conversion

library(tidyverse)
source("full-tasr-isd/cleaning_functions.R")

file_names <- list.files("full-tasr-isd/annual_length_data/length_2015-2018",
                         pattern = "*.csv",
                         full.names = TRUE)

column_names <- list()
for(i in 1:length(file_names)){
  cnames <- names(read.csv(file_names[i]))
  column_names[[i]] <- cnames
}
column_names

target_names_all <- c("Stream",  "Fraction", 
                      "L1", "L2",
                      "L3", "L4", "L5", "L6", "L7", 
                      "L8", "L9", "L10", "L12", "L13",
                      "L14", "L15", "L16", "L17", "L18",
                      "L19", "L20", "Subsample", "Abundance")

target_names_15_16 <- c("Site", "Order", "Family", "Genus", "X..of.Surbers")
target_names_17_18 <- c("Rep", "Taxa", "NoSurbers")

for(i in seq_along(column_names)){
  print(setdiff(target_names_all, column_names[[i]]))
  #print(setdiff(column_names[[i]], target_names))
}

# 2015 and 2016
# filter out Abundance == 0 or NA
# Read in order, family, genus
# read in data
dat_list <- list()
for(i in 1:2){
  temp_dat <- read.csv(file_names[i])
  temp_dat <- temp_dat |>
    filter(Abundance > 0,
           !is.na(Abundance))
  temp_dat <- temp_dat[,names(temp_dat) %in% c(target_names_all, target_names_15_16)]
  temp_dat <- temp_dat |>
    mutate(Genus = case_when(Genus == "Unk" ~ Family,
                             .default = Genus),
           Genus = case_when(Genus == "Unk" ~ Order,
                             .default = Genus),
    )
  temp_dat$Taxa <- temp_dat$Genus
  temp_dat$NoSurbers <- temp_dat$X..of.Surbers
  temp_dat$Rep <- 1
  temp_dat$Year <- substr(
    list.files("full-tasr-isd/annual_length_data/length_2015-2018",
               pattern = "*.csv")[i], 1, 4)
  temp_dat <- temp_dat[,names(temp_dat) %in% c("Year", "Site", target_names_17_18, target_names_all)]
  temp_dat <- temp_dat[,c("Year", "Site", target_names_17_18, target_names_all)]
  dat_list[[i]] <- temp_dat
}
dat_15_16 <- bind_rows(dat_list)
names(dat_15_16)
dat_15_16 |> distinct(NoSurbers)
dat_15_16 |> distinct(Subsample)
dat_15_16 |> distinct(NoSurbers, Subsample)

dat_list <- list()
for(i in 3:4){
  temp_dat <- read.csv(file_names[i])
  temp_dat <- temp_dat[,names(temp_dat) %in% c(target_names_all, target_names_17_18)]
  temp_dat$Year <- substr(
    list.files("full-tasr-isd/annual_length_data/length_2015-2018",
               pattern = "*.csv")[i], 1, 4)
  temp_dat$Site <- NA
  temp_dat <- temp_dat[,names(temp_dat) %in% c("Year", "Site", target_names_17_18, target_names_all)]
  temp_dat <- temp_dat[,c("Year", "Site", target_names_17_18, target_names_all)]
  dat_list[[i]] <- temp_dat
}
dat_17_18 <- bind_rows(dat_list)
dat <- bind_rows(dat_15_16, dat_17_18)
as_tibble(dat)
# need to deal with subsamples or unequal sample size?
dat |> distinct(Subsample)


# stream names
distinct(dat, Stream, Site, Year) |> arrange(Stream, Year)

# Taxa names
distinct(dat, Taxa) |> arrange(Taxa)

dim(dat)
dat_cleaned <- clean_2015_2018(dat)
dim(dat_cleaned)

# Taxa names
distinct(dat_cleaned, Taxa) |> arrange(Taxa)

# Duplicated rows ####
check_duplicates(dat_cleaned) 
# I think this is mostly due to the upper/lower streams
check_duplicates(dat_cleaned) |>
  ungroup() |>
  distinct(Year, Stream)
# removing these rows for now
# Talk to Lusha and see how to proceed

remove_target <- check_duplicates(dat_cleaned)
dim(dat_cleaned)
dim(anti_join(dat_cleaned, remove_target))
dat_cleaned <- anti_join(dat_cleaned, remove_target)

# pivot data and sample "missing" measurements
dim(dat_cleaned)
# pivot_ function also removes rows where Abundance is missing
# Ask Lush about missing totals
set.seed(8789)
dat_sim <- pivot_and_sample(dat_cleaned)
dim(dat_sim)
dat_sim
l_dat <- dat_for_lw(dat_sim)
l_dat
l_dat |> distinct(Taxa) |> arrange(Taxa) |> pull()

l_dat |> 
  clean_taxa_names() |>
  distinct(Taxa) |>
  arrange(Taxa) |>
  pull()

l_dat <- l_dat |>
  clean_taxa_names()

saveRDS(l_dat, "full-tasr-isd/cleaned_length_data/2015-2018.rds")
