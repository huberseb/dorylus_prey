library(ggplot2)
library(ggthemes)
library(patchwork)

#### Biometrics ####

# read in file
# dec = comma!
all_df <- read.delim(
  file   = "all_ants_all_prey.csv",
  sep  = ",",
  dec = "." # <-- check if eng with "point "." or german with "," 
)


# Plot ant mandible vs. body length for BOTH species
ggplot(all_df, aes(x = ant_length, y = ant_width, color = ant_species)) +
  geom_point(alpha = 0.5, size = 1) +             # slightly transparent points
  geom_smooth(method = "lm", se = TRUE, linewidth = 1) +  # linear regression line 
  theme_clean(base_size = 14) +
  labs(
    x = "Ant body length (mm)",
    y = "Ant body width (mm)",
    color = "Species",
    title = "Relationship between body length and width in Dorylus ants"
  )

# Fuck around - Want to change coulours but dosnt work for now...
ggplot(all_df, aes(x = ant_length, y = ant_width, color = ant_species)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.2) +
  scale_color_manual(values = c("D. wilverthi" = "firebrick3",
                                "D. sjostedti" = "steelblue3")) +
  theme_clean(base_size = 14) +
  labs(
    x = "Ant body length (mm)",
    y = "Ant body width (mm)",
    color = "Species",
    title = "Body length–width relationship in two Dorylus species"
  )

# Plot ant mandible vs. body length for EACH species
ggplot(all_df, aes(x = ant_length, y = ant_width)) +
  geom_point(color = "darkred", alpha = 0.7, size = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ ant_species) +
  theme_clean(base_size = 14) +
  labs(
    x = "Ant body length (mm)",
    y = "Ant body width (mm)",
    title = "Length–width relationship by Dorylus species"
  )
