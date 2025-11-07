# we will plot stuff here


library(ggplot2)
library(ggthemes)

#### Biometrics ####

# read in file
# dec = comma!
raw_multi <- read.delim(
  file   = "multiple_carriers.csv",
  sep  = ";",
  dec = "," # <-- very important for comma decimals
)


# Plot ant mandible vs. body length 


ggplot(df, aes(x = length, y = width, color = ant.species)) +
  geom_point(alpha = 0.7, size = 2) +             # Punkte mit leichter Transparenz
  geom_smooth(method = "lm", se = TRUE, linewidth = 1) +  # lineare Regressionslinie mit Konfidenzintervall
  theme_clean(base_size = 14) +
  labs(
    x = "Ant body length (mm)",
    y = "Ant body width (mm)",
    color = "Species",
    title = "Relationship between body length and width in Dorylus ants"
  )