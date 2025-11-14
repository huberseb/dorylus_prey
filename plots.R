library(ggplot2)
library(ggthemes)
library(patchwork)

#=============================================================================#
                            #### Biometrics #### 
#=============================================================================#

# read in file
# dec = comma!
all_df <- read.delim(
  file   = "all_ants_all_prey.csv",
  sep  = ",",
  dec = "." # <-- check if eng with "point "." or german with "," 
)


# Ant length vs. ants width for BOTH species###########################
ant_length_width <- ggplot(all_df, aes(x = ant_length, y = ant_width, color = ant_species)) +
  geom_point(alpha = 0.3, size = 1.5) +             # slightly transparent points
  geom_smooth(aes(color = NULL),    # remove species grouping for the line
              method = "lm", 
              se = TRUE, 
              linewidth = 1,
              color = "grey9") + 
  scale_x_log10() +     
  scale_y_log10() + 
  scale_colour_manual(values = c("D. wilverthi" = "#56B4E9",
                                 "D. sjostedti" = "#E69F00")) +
  theme_clean(base_size = 14) +
  theme(legend.position = "bottom", plot.title = element_text(size = 14)) +
  labs(
    x = "Ant body length (mm)",
    y = "Ant body width (mm)",
    color = "Species",
    title = "Relationship between body length and width in Dorylus ants"
  )

jpeg(filename = "Relationships length width ants.jpg", width = 20, 
     height = 18, units = "cm", res = 300)
ant_length_width
dev.off()

#Ant length for both species##################################################
ggplot(all_df, aes(x = ant_length, colour = ant_species, fill =  ant_species)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  scale_fill_manual(values = c("D. wilverthi" = "#56B4E9",
                                "D. sjostedti" = "#E69F00")) +
  scale_colour_manual(values = c("D. wilverthi" = "#56B4E9",
                               "D. sjostedti" = "#E69F00")) +
  guides(colour = "none") +  
  theme_clean(base_size = 14) +
  theme(plot.title = element_text(hjust = 1)) + 
  labs(
    x = "Ant body length [mm]",
    y = "Frequency",
    fill = "Dorylus species",
    title = "Body length in two Dorylus species"
  )

#Ant widht vs. body length for EACH species###########################
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

hist(df$prey_weight)


#=============================================================================#
                      #### Linear Regression models #### 
#=============================================================================#



#Ant weight vs prey weight (single workers)###################################
#This is the plot showing the lm3
#Regression line shoulde be for bothe species sinces species has no sig. effect

ggplot(df_single, aes(x = log_ant_weight, y = log_prey_weight,
                      color = ant_species)) +
#sets points  
  geom_point(alpha = 0.7, size = 0.7) +  
  
#regression line, se sets confidence area 
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7) +

#add horizontel reference line   
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, color = "red") +
  
#set colour scheme 
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +

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

  
#Relative load vs. Carrier size (single workers)###############################

  ggplot(df_single, aes(x = log10(ant_weight), y = log10(relativ_loading),
                        color = ant_species)) +
    #sets points  
    geom_point(alpha = 0.7, size = 0.7) +  
    
    #regression line, se sets confidence area 
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
    
    #add horizontel reference line   
    #geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, color = "red") +
    
    #set colour scheme 
    scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  
  #Add description and legend
  labs(
    x = "Ant weight",
    y = "Relative loading",
    color = "Dorylus"
  ) +
    theme_classic() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12)
    )
