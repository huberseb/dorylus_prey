library(ggplot2)
library(effectsize)
library(ggeffects)
library(tidyvers)
library(tidyverse)
library(dplyr)
library(lmtest)
library(performance)
library(lme4) #GLMM
library(DHARMa)

#-----------------------------------------------------------------------------#
                    ####1.Prey Spectra by Ant Species ####
#-----------------------------------------------------------------------------#


#-----------------------------------------------------------------------------#
            ####2.Dismembered - Preparation for transport####
#-----------------------------------------------------------------------------#
##### 2.1 Clean the dismembered column #####
#reading in data
df <- read.delim(
  file   = "prey_spectra_final.csv",
  sep  = ";",
  dec = "," # <-- very important for comma decimals
)

# Normalize formatting (lowercase + trim white space)
df$dismembered <- tolower(trimws(df$dismembered))

# Keep only valid states ("yes" / "no")
df_dis <- df %>%
  filter(dismembered %in% c("yes", "no"))

table(df_dis$dismembered)

# Make sure its a factor
df_dis$dismembered <- factor(df_dis$dismembered, levels = c("no", "yes"))
df_dis$ant_species <- factor(df_dis$ant_species)
df_dis$forest_type <- factor(df_dis$forest_type)

# Add important factors 
df_dis["prey_shape"] <- df_dis$prey_length / df_dis$prey_width
df_dis["prey_area"] <- df_dis$prey_length * df_dis$prey_width


##### 2.2 Dismembered stat vs. Prey Shape (prey_length / prey width)####


# Boxplot of prey_shape by dismembered, faceted by species
# Just to get a feeling
ggplot(df_dis, aes(x = dismembered,
                   y = log10(prey_shape))) +
  geom_boxplot() +
  facet_wrap(~ ant_species) +
  labs(
    x = "Dismembered (no / yes)",
    y = "log10(prey_length / prey_width)",
    title = "Prey shape vs. dismembering state by ant species"
  )

# GLM with binomial factors
# Logistic regression dismembered ~ prey shape * species + forest_type + prey_order
mod_shape_full <- glm(
  dismembered ~ log10(prey_shape) * ant_species + forest_type + prey_order,
  data   = df_dis,
  family = binomial
)

summary(mod_shape_full)
anova(mod_shape_full, test = "Chisq") #Likelihood-ratio tests for model terms



##### 2.2 Dismembered stat vs. Prey Area (prey_length * prey width)####

# Boxplot of prey_area by dismembered, faceted by species
# Just to get a feeling
ggplot(df_dis, aes(x = dismembered,
                   y = log10(prey_area))) +
  geom_boxplot() +
  facet_wrap(~ ant_species) +
  labs(
    x = "Dismembered",
    y = "log10(prey_area)",
    title = "Prey Area vs. dismembering state by ant species"
  )

# GLM with binomial factors
# Logistic regression dismembered ~ prey shape * species + forest_type + prey_order
mod_area_full <- glm(
  dismembered ~ log10(prey_area) * ant_species + forest_type + prey_order,
  data   = df_dis,
  family = binomial
)

summary(mod_area_full)
anova(mod_area_full, test = "Chisq") #Likelihood-ratio tests for model terms


##### 2.3 Threshold of dismembering by a GLMM #####
library(lme4)

#logit-GLMM with log-area, species, Raid as Random Intercept
m_area_int <- glmer(
  dismembered ~ log10(prey_length * prey_width) * ant_species +
    (1 | raid_ID),
  data   = df_dis,
  family = binomial(link = "logit")
)

summary(m_area_int)

#optional Model with added Random Slope
m_area_slope <- glmer(
  dismembered ~ log10(prey_length * prey_width) * ant_species +
    (1 + log10(prey_length * prey_width) | raid_ID),
  data   = df_dis,
  family = binomial(link = "logit")
)

summary(m_area_slope) # <- its the fitter model, AIC and BIC are smaller

#used m_area_slope for further calculations
plot(m_area_slope)
testDispersion(m_area_slope)
testUniformity(m_area_slope)

#Calculating 50% Treashold with tested model form above
#since the summary showed no differences between species we calculate a combined

fe <- fixef(m_area_slope) #gives out fixed effects of the model 

b0 <- fe["(Intercept)"]
b1 <- fe["log10(prey_length * prey_width)"]

# log10(area) at 50% probability
x_thresh_comb <- -b0 / b1

# taking the area out (mm^2)
area_thresh_comb <- 10^x_thresh_comb

x_thresh_comb
area_thresh_comb # <- this is our threshold/result


# ---------------------------------------------------------------------------#
                        #generate predictions from GLMM
#----------------------------------------------------------------------------#
pred <- ggpredict(
  m_area_slope,
  terms = c("log10(prey_length * prey_width)", "ant_speciesD. wilverthi")
)

                              # 50%-Threshold 
fe <- fixef(m_area_slope)
b0 <- fe["(Intercept)"]
b1 <- fe["log10(prey_length * prey_width)"]
x_thresh_comb <- -b0 / b1  # log10(area)
area_thresh_comb <- 10^x_thresh_comb  # area in mm^2


#Plot GLMM <- NOT WORKING YET

ggplot(pred, aes(x = x, y = predicted, colour = group)) +
  geom_line(size = 1.1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group),
              alpha = 0.2, colour = NA) +
  geom_point(data = df_dis,
             aes(x = log10(prey_length * prey_width),
                 y = dismembered),
             inherit.aes = FALSE,
             alpha = 0.15, size = 1, colour = "grey30") +
  geom_vline(xintercept = x_thresh_comb,
             linetype = "dashed",
             colour = "red",
             size = 1) +
  labs(
    x = "log10(Prey area [mm²])",
    y = "Predicted probability of dismemberment",
    colour = "Species",
    fill = "Species",
    title = "Size-dependent dismemberment probability in Dorylus ants",
    subtitle = paste0("50% threshold at ~", round(area_thresh_comb, 2), " mm²")
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )



#-----------------------------------------------------------------------------#
                          ####3.Life stage of Pey####
#-----------------------------------------------------------------------------#