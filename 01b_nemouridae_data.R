# Simulating Nemouridae data from subsample of measured individuals

library(tidyverse)

upperDir <- "Nemouridae-data" 

# vector of all of the paths to csv's within folder
meas <- list.files(upperDir, full.names=TRUE) 
# just the file name of each csv
filename <- list.files(upperDir) 

# nem_filename <- list.files("Nemouridae-data") 
# raw_filename <- list.files("raw_dat") 
# nem_filename
# setdiff(nem_filename, raw_filename)
#setdiff(raw_filename, nem_filename)

# number of files?
length(filename)

#### Loop Reading in csv's ####
dataset <- NULL #Empty dataset


start_time <- Sys.time()

for (i in 1:length(meas)){ 
  #length(meas) = number of total files within folder
  
  # if errors occur when reading in files, you may need to add an empty line to the csv files. 
  # Uncomment the following once 
  
  # adding empty line to csv files
  #cat("\n", file = meas[i], append = TRUE)
  
  
  # creating temporary dataset for file
  temp_dataset <- read.csv(meas[i], header = TRUE) 
  temp_dataset <- temp_dataset[1:4] # forcing file to only have 4 columns. When measurements are exported into excel, Adobe sometimes adds an additional column
  temp_dataset$file <- filename[i] # adding new colum "file" to temp data set
  dataset <- rbind(dataset,
                   temp_dataset) #row binding from temp dataset to full dataset
  rm(temp_dataset) #clearing out temp data for next file
}

# time for code to run?
end_time <- Sys.time()
end_time - start_time
dataset

# checking that data matches
dataset |>
  group_by(file) |>
  count()

dataset |>
  filter(file == "DELT_01_LG_2018_08_02_382.csv") |>
  pull(Value)

#  total counts for sub-samples -------------------------------------------
# This information came from an email from Lusha on 2026-04-19
# files that were subsampled
sub_files <- c("DELT_01_LG_2018_08_02_382.csv",
               "DELT_01_SM_2018_08_02_383.csv",
               "DELT_02_LG_2018_08_02_384.csv",
               "DELT_02_SM_2018_08_02_385.csv")


tot_df <- tibble(
  file = sub_files,
  Label = c("Zapada",
            "Zapada",
            "Zapada",
            "Zapada"),
  tot = c(117,
          103,
          31,
          107)
)

# join the full data set to the total counts
tot_df_meas <- left_join(tot_df, dataset)

tot_df_meas |>
  group_by(file) |>
  count()

# re-sample the individual measurements ####
nemouridae_sim <- tot_df_meas |>
  group_by(file) |>
  mutate(n = n(), 
         w = round(tot / n)) |>
  #ungroup() |>
  rowwise() |>
  sample_n(w, 
           replace = T)


# combining "full" and simulated nemouridae data --------------------------
nemouridae_emp <- dataset |>
  filter(!file %in% sub_files)
# check that the right files were removed/kept
nemouridae_emp |>
  group_by(file) |>
  count()

# combine empirical and sub-sampled data
nemouridae_dat <- bind_rows(nemouridae_emp,
                            nemouridae_sim) |>
  as_tibble()

# separate file column ####
nemouridae_dat <- nemouridae_dat |> 
  as_tibble() |>
  select(-tot, -n, -w) |>
  separate(file, c("site", "rep", "size_fraction", 
                   "year", "month", "day", 
                   "photo_number"),
           sep = "_", remove = FALSE)
nemouridae_dat

# fix day-range
nemouridae_dat <- nemouridae_dat |> 
  separate(day, 
           c("day", "extra_day"),
           sep = "-", 
           remove = TRUE) |>
  select(-extra_day)

# format date ####
nemouridae_dat <- nemouridae_dat %>%
  mutate(
    year = as.Date(year, format = "%Y"),
    year = year(year),
    month = month(as.numeric(month)),
    day = as.Date(day, format = "%d"),
    day = day(day))
nemouridae_dat

# change label ####
nemouridae_dat <- nemouridae_dat |>
  mutate(Label = "Nemouridae") 

# add lw coefs ####
# read in csv with length weight equations
lw_coef <- read.csv("LW_coeffs.csv")

lw_cols <- c("taxon",
             "a",
             "b",
             "L_units",
             "dw_units",
             "formula_type",
             "formula")

# check formulas in LW for taxa in arkansas
equations <- lw_coef[lw_coef$taxon %in% unique(nemouridae_dat$Label),]


# join the nemouridae data with lw coef values ####
nem_eq <- merge(nemouridae_dat, lw_coef[,lw_cols], 
                 by.x = "Label",
                 by.y = "taxon", 
                 all.x = FALSE)


# estimate dw based on formula type ####
nem_eq <- nem_eq %>% 
  mutate(dw = case_when(
    formula_type == 1 ~ a * Value^b,
    formula_type == 2 ~ exp(a + b * log(Value))))


as_tibble(nem_eq)
# make nemouridae match full data ####
# remove any dry weight values < or = to 0
nem_dw <- nem_eq %>%
  filter(dw >0) %>%
  select(site, rep, year, month, file, Label, Value, dw) |>
  as_tibble()

nem_dw
# save nemouridae data
write_rds(nem_dw, "derived_data/TASR_nemouridae_dw.rds")
