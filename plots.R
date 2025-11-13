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


# Plot ant mandible vs. body length for BOTH species###########################
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

# Plot ant mandible vs. body length for EACH species###########################
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

#=============================================================================#
                      #### Linear Regression models #### 
#=============================================================================#

#Ant weight vs prey weight (single workers)###################################
#This is the plot showing the lm1
#Regression line shoulde be for bothe species sinces species has no sig. effect

ggplot(df_single, aes(x = log_ant_weight, y = log_prey_weight,
                      color = ant_species)) +
#sets points  
  geom_point(alpha = 0.7, size = 0.7) +  
  
#regression line, se sets confidence area 
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +

#add horizontel reference line   
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, color = "red") +
  
#set colour scheme 
  scale_color_manual(values = c("#E69F00", "#56B4E9")) 

#Add description and legend
  labs(
    x = "Ant weight (log10, mg)",
    y = "Prey weight (log10, mg)",
    color = "Dorylus"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )
  #TO DO: CONVERT g TO mg FOR NICER SCALING. (TILL THEN LEGEND WRONG)
  #       Check for used model (is it the right one?)

