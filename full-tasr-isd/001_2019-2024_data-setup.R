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

file_names <- list.files("full-tasr-isd/annual_length_data",
                         pattern = "*.csv",
                         full.names = TRUE)

column_names <- list()
for(i in 1:length(file_names)){
  cnames <- names(read.csv(file_names[i]))
  column_names[[i]] <- cnames
}
column_names

target_names <- c("Year", "Stream", "Rep",
                  "Fraction", "Taxa", "L1", "L2",
                  "L3", "L4", "L5", "L6", "L7", 
                  "L8", "L9", "L10", "L12", "L13",
                  "L14", "L15", "L16", "L17", "L18",
                  "L19", "L20", "Subsample", "Abundance", "Stage")

for(i in seq_along(column_names)){
  print(setdiff(target_names, column_names[[i]]))
  print(setdiff(column_names[[i]], target_names))
}

# read in data
dat_list <- list()
for(i in seq_along(file_names)){
  temp_dat <- read.csv(file_names[i])
  temp_dat <- temp_dat[,names(temp_dat) %in% target_names]
  temp_dat <- temp_dat[,target_names]
  # fixing "Metabarcoding" rep in two streams
  # the samples were switched out, but the metabarcoding sample was collected and processed in the same way. 
  if(temp_dat$Year[1] == 2024) {
    temp_dat <- temp_dat |>
      mutate(
        Rep = case_when(
          Stream == "Grizzly" & Rep == "Metabarcoding" ~ "1",
          Stream == "NF Teton Creek" & Rep == "Metabarcoding" ~"2",
          .default = Rep),
        Rep = as.numeric(Rep))
  }
  dat_list[[i]] <- temp_dat
}
dat <- bind_rows(dat_list)

# need to deal with subsamples?
dat |> distinct(Subsample)
# no, all fully picked (Subsample = 1) or NA - I assume that means fully picked?

# stream names
distinct(dat, Stream, Year) |> arrange(Stream, Year)

# Taxa names
distinct(dat, Taxa) |> arrange(Taxa)

dim(dat)
dat_cleaned <- clean_2019_2024(dat)
dim(dat_cleaned)

# Taxa names
distinct(dat_cleaned, Taxa) |> arrange(Taxa)

check_duplicates(dat_cleaned) 
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

saveRDS(l_dat, "full-tasr-isd/cleaned_length_data/2019-2024.rds")
