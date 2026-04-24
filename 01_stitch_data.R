library(tidyverse)

# Download all files from one drive intoa local folder
# 
upperDir <- "raw_dat" 

# vector of all of the paths to csv's within folder
meas <- list.files(upperDir, full.names=TRUE) 
# just the file name of each csv
filename <- list.files(upperDir) 

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



# check data for errors
dataset
names(dataset)
head(dataset)

# measurement units ####
# all units
dataset %>%
  distinct(Unit)


dataset %>%
  filter(
         Unit == "in") %>% 
  group_by(file) %>%
  count() %>%
  ungroup() %>%
  distinct(n)
  

head(dataset)

dataset <- dataset %>% 
  as_tibble() %>%
  separate(file, c("site", "rep", "size_fraction", 
                   "year", "month", "day", 
                   "photo_number"),
           sep = "_", remove = FALSE)

### Change days to one number see if that fixes month prob below
### it does not :(

distinct(dataset, day)
distinct(dataset, month)
as_tibble(dataset)


#View(dataset)
dataset <- dataset %>%
  mutate(
    year = as.Date(year, format = "%Y"),
    year = year(year),
    month = month(as.numeric(month)),
    day = as.Date(day, format = "%d"),
    day = day(day)) |>
  as_tibble()
dataset

# distinct(dataset,
#          site, 
#          date)
# dataset %>% distinct(year, site, rep) %>% arrange(year, site)
# label typos ####
sort(unique(dataset$Label))
# lots of typos to fix

dataset %>%
  filter(Label == "") %>%
  distinct(file)

dataset %>%
  filter(Label == "scale") %>%
  distinct(file)


# fix Label Typos ####

sort(unique(dataset$Label))

# finish fixing typos!!
# save as dat_clean once cleaned
dat_clean <- dataset %>%
  mutate(Label = case_when(
    ##++++++++++++++++####
    ##+ Change lednia and zapada to "Nemouridae"
    ##+ 
    Label == "simullidae" ~ "Simuliidae",
    Label == "baetidae" ~ "Baetidae",
    Label == "beetle?" ~ "Elmidae",
    Label == "chirnomidae" ~ "Chironomidae",
    Label == "chironomida" ~ "Chironomidae",
    Label == "chironomidae" ~ "Chironomidae",
    Label == "chironomidaea" ~ "Chironomidae",
    Label == "chloroperlidae" ~ "Chironomidae",
    Label == "chrionomidae" ~ "Chironomidae",
    Label == "Chrionomidae" ~ "Chironomidae",
    Label == "chrironomidae" ~ "Chironomidae",
    Label == "chronomidae" ~ "Chironomidae",
    Label == "limnephilidae" ~ "Limnephilidae",
    Label == "oligiochaeta" ~ "Oligochaeta",
    Label == "oligo" ~ "Oligochaeta",
    Label == "oligochaeta" ~ "Oligochaeta",
    Label == "oloigochaeta" ~ "Oligochaeta",
    Label == "perlodidae" ~ "Perlodidae",
    Label == "perlolidae" ~ "Perlodidae",
    Label == "planaria" ~ "Planaria",
    Label == "rhyacophilia" ~ "Rhyacophilidae",
    Label == "rhyacophilidae" ~ "Rhyacophilidae",
    Label == "simuiidae" ~ "Simuliidae",
    Label == "simuliidae" ~ "Simuliidae",
    Label == "simullidae" ~ "Simuliidae",
    Label == "tipulidae" ~ "Tipulidae",
    Label == "pt 2" ~ "Chironomidae",
    Label == "plecoptera" ~ "Chloroperlidae",
    .default = Label
  )) 

dat_clean %>%
  filter(Label == "Non-insects")


# read in csv with length weight equations
lw_coef <- read.csv("LW_coeffs.csv")

# which taxa don't have lw coefs?
no_coef <- setdiff(unique(dat_clean$Label), unique(lw_coef$taxon))

# should be nothing after fixing typos
sort(no_coef)

# Using surrogate equations for these taxa
# Morley et al. 2020 Plos One
# Psychodidae --> Empididae
# Uenoidae --> Limnephilidae
# Hydroptilidae --> Rhyacophilidae
# dat_clean <- dat_clean %>%
#   mutate(Label = case_when(
#     Label == "Psychodidae" ~ "Empididae",
#     Label == "Uenoidae" ~ "Limnephilidae",
#     Label == "Hydroptilidae" ~"Rhyacophilidae",
#     .default = Label
#   )) 

no_coef2 <- setdiff(unique(dat_clean$Label), unique(lw_coef$taxon))

sort(no_coef2)

# percent of data that does not have lw coeffs?
(nrow(dat_clean[dat_clean$Label %in% no_coef2,])/ nrow(dat_clean))*100

# estimate dry weight (dw) ####

# make string vector of column names we care about
lw_cols <- c("taxon",
             "a",
             "b",
             "L_units",
             "dw_units",
             "formula_type",
             "formula")

# check formulas in LW for taxa in arkansas
equations <- lw_coef[lw_coef$taxon %in% unique(dat_clean$Label),]


# join the ark data with lw coef values
fulldat <- merge(dat_clean, lw_coef[,lw_cols], 
                 by.x = "Label",
                 by.y = "taxon", 
                 all.x = FALSE)

# dimensions of the ark and merged data
dim(dataset)
dim(fulldat)

# different lw formulas
unique(fulldat$formula_type)

# estimate dw based on formula type
fulldat <- fulldat %>% 
  mutate(dw = case_when(
    formula_type == 1 ~ a * Value^b,
    formula_type == 2 ~ exp(a + b * log(Value))))


names(fulldat)
# remove any dry weight values < or = to 0
dw <- fulldat %>%
  filter(dw >0) %>%
  select(site, rep, year, month, file, dw)

nrow(fulldat)
nrow(dw)

names(dw)

as_tibble(dw)


# save data with estimated dry weights
write_rds(dw, "derived_data/TASR_dw.rds")
# change file path name^^^
# make a derived data folder to save in








