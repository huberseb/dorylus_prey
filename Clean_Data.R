library(readr)
library(dplyr)
library(tidyr)
library(stringr)

#### read data ####
# read multiple carriers csv file
# dec = comma!
raw_multi <- read.delim(
  file   = "multiple_carriers.csv",
  sep  = ";",
  dec = "," # <-- very important for comma decimals
)

# read prey spectra csv file
# dec = comma!
raw_prey <- read.delim(
  file   = "prey_spectra_all_carriers.csv",
  sep  = ";",
  dec = "," # <-- very important for comma decimals
)


#### Make changes to "multiple" dataframes ####

# make long format; one row = one ant
long_multi <- pivot_longer(
  data = raw_multi,
  cols = matches("^ant_(length|width|weight)\\d+$"),
  names_to   = c(".value","ant_index"),
  names_pattern = "^ant_(length|width|weight)(\\d+)$"
) 

# remove all rows where NA is present (max 14 ants per tube --> many rows with NA)
long_multi_clean <- drop_na(long_multi)

# rename to sum to not confuse with other tables later
# do this just to revert it again in 5 steps....
raw_prey <- rename(raw_prey, sum_ant_length = ant_length)
raw_prey <- rename(raw_prey, sum_ant_width = ant_width)
raw_prey <- rename(raw_prey, sum_ant_weight = ant_weight)

# merge prey and multi to add prey info to multiple carrier ant info
df_merge_multi <- merge(raw_prey, long_multi_clean, by="tube_ID")

# take from df_merge all columns except...
# multiple prey info from both dataframes caused duplicates
df_merge_multi <- df_merge_multi %>% select(!c(prey_length.y, prey_weight.y, prey_width.y))
# rename remainders; to have same column name as future dfs
df_merge_multi <- rename(df_merge_multi, prey_width = prey_width.x)
df_merge_multi <- rename(df_merge_multi, prey_length = prey_length.x)
df_merge_multi <- rename(df_merge_multi, prey_weight = prey_weight.x)
df_merge_multi <- rename(df_merge_multi, ant_number = ant_number.x)

# from raw_prey, only keep those that were singles
raw_prey_singles <- raw_prey %>% filter(ant_number == 1)

# now, we revert the steps we took before. might be possible to not do any of this or the ones before
# TODO: check if possible not to rename to sum_
raw_prey_singles$length = raw_prey_singles$sum_ant_length
raw_prey_singles$width = raw_prey_singles$sum_ant_width
raw_prey_singles$weight = raw_prey_singles$sum_ant_weight

# set ant_index for singles to 1 (obacht: make it a character for next step)
raw_prey_singles$ant_index <- 1
raw_prey_singles <- raw_prey_singles %>% mutate(ant_index = as.character(ant_index))

# combine singles and multiples
combined_df <- bind_rows(raw_prey_singles, df_merge_multi)

#loose ant_number.y
# TODO: check if neccesary, when not called ant_number.y, but ant_number
combined_df <- combined_df %>% select(!c(ant_number.y))


# make ant_index numeric, so we can sort
combined_df <- combined_df %>% mutate(ant_index = as.numeric(ant_index))

####Final clean up ####

# sort colums to a more fitting order 25
combined_df_sorted <- combined_df[, c("tube_ID", "ant_index", "length", "width", "weight", "ant_number",
 "sum_ant_length", "sum_ant_width", "sum_ant_weight", "ant_species",
"prey", "prey_length", "prey_width", "prey_weight", "dismembered.", "life_stage", "holometabola",
 "prey_phylum", "prey_class", "prey_order",
"prey_family", "prey_genus", "raid_ID", "colony_ID", "forest_type")]

# name all colums fitting 
combined_df_sorted <- rename(combined_df_sorted, ant_width = width)
combined_df_sorted <- rename(combined_df_sorted, ant_length = length)
combined_df_sorted <- rename(combined_df_sorted, ant_weight = weight)
combined_df_sorted <- rename(combined_df_sorted, dismembered = dismembered.)
combined_df_sorted <- rename(combined_df_sorted, total_ant_number = ant_number)


##### save final dataframe as .csv #####
# think of nice filename once issues are fixed
write.csv(combined_df_sorted, "all_ants_all_prey.csv")


