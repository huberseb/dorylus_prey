library(readr)
library(dplyr)
library(tidyr)
library(stringr)

# read multiple carriers csv file
# dec = comma!
raw_multi <- read.delim(
  file   = "multiple_carriers.csv",
  sep  = ";",
  dec = "," # <-- very important for comma decimals
)

# make long format; one row = one ant
long <- pivot_longer(
  data = raw,
  cols = matches("^ant_(length|width|weight)\\d+$"),
  names_to   = c(".value","ant_index"),
  names_pattern = "^ant_(length|width|weight)(\\d+)$"
) 

# remove all rows where NA is present (max 14 ants per tube --> many rows with NA)
long_clean_multi <- drop_na(long)


# read prey spectra csv file
# dec = comma!
raw_prey <- read.delim(
  file   = "prey_spectra_all_carriers.csv",
  sep  = ";",
  dec = "," # <-- very important for comma decimals
)

raw_prey <- rename(raw_prey, sum_ant_length = ant_length)
raw_prey <- rename(raw_prey, sum_ant_width = ant_width)
raw_prey <- rename(raw_prey, sum_ant_weight = ant_weight)

df_merge <- merge(raw_prey, long_clean_multi, by="tube_ID")

df_merge_long_first <- merge(long_clean_multi,raw_prey, by="tube_ID")

raw_prey_singles <- raw_prey %>% filter(ant_number == 1)

combined_df <- bind_rows(raw_prey_singles, df_merge)
