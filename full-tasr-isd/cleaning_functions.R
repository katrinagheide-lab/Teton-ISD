clean_2019_2024 <- function(df){
  df <- df |>
    filter(!Stage %in% c("adult", "pupae"),
           !str_detect(Taxa, regex("adult", ignore_case = TRUE)),
           !str_detect(Taxa, regex("pupa", ignore_case = TRUE)))
  return(df)
}

clean_2015_2018 <- function(df){
  df <- df |>
    filter(
           !str_detect(Taxa, regex("adult", ignore_case = TRUE)),
           !str_detect(Taxa, regex("pupa", ignore_case = TRUE)))
  return(df)
}

check_duplicates <- function(df){
  df |>
    group_by(Stream, Year, Rep, Taxa, Fraction) |>
    count() |>
    filter(n > 1)
}

pivot_and_sample <- function(df){
  if(any(is.na(df$Abundance))){
    df <- df|>
      filter(!is.na(Abundance))
    message("NA Abundances removed")
  }
  out <- df |>
    pivot_longer(L1:L20) |>
    filter(!is.na(value)) |>
    select(Stream, Year, Taxa, Abundance, Fraction, Rep, name, value) |>
    group_by(Stream, Year, Taxa, Rep, Fraction) |>
    mutate(n = n(),
           replace = Abundance>n,
           prop = Abundance / n) |>
    group_modify(~ slice_sample(.x,
                                prop = unique(.x$prop),
                                replace = unique(.x$replace)))
  return(out)
}

dat_for_lw <- function(df){
  out <- df |>
    ungroup() |>
    select(Year, Stream, Taxa, l_mm = value)
  return(out)
}

# clean Taxa names ####
clean_taxa_names <- function(df){
  out <- df |>
    mutate(Taxa = ifelse(str_detect(Taxa, "Rithrogena"), "Rhithrogena", Taxa),
           Taxa = ifelse(str_detect(Taxa, "Megarcys"), "Megarcys", Taxa),
           Taxa = ifelse(str_detect(Taxa, "Baetis"), "Baetis", Taxa),
           Taxa = ifelse(str_detect(Taxa, "Homphylax"), "Homophylax", Taxa),
           Taxa = ifelse(str_detect(Taxa, "Rhyacophila"), "Rhyacophila", Taxa),
           Taxa = ifelse(str_detect(Taxa, "Molophilus did"), "Limoniidae", Taxa),
           Taxa = ifelse(str_detect(Taxa, "Non-Tany"), "Chironomidae", Taxa),
           Taxa = ifelse(str_detect(Taxa, "Simuliidae"), "Simuliidae", Taxa),
           Taxa = ifelse(str_detect(Taxa, "Arctotipula"), "Tipula", Taxa),
           Taxa = ifelse(str_detect(Taxa, "Tricoptera"), "Trichoptera", Taxa),
           Taxa = ifelse(str_detect(Taxa, "diptera"), "Diptera", Taxa))
  return(out)
}
