library(ggplot2)
library(ggthemes)
library(patchwork)
library(effectsize)
library(tidyverse)
library(dplyr)
library(lmtest)
library(performance) #e.g. for ICC calculations
library(lme4) #for LMM´s
library(DHARMa)

#reading in data
df <- read.delim(
  file   = "all_ants_all_prey.csv",
  sep  = ",",
  dec = "." # <-- very important for comma decimals
)
df$ant_species <- factor(df$ant_species)
df$forest_type <- factor(df$forest_type)

#==============================================================================#
                      ####  Preparation fo Data ####
#==============================================================================#
#####  Testing Distribution ##### 
#==============================================================================#
#Shapiro-Wilk test 
shapiro.test(df$ant_length)
shapiro.test(df$ant_width)
shapiro.test(df$ant_weight)
shapiro.test(df$prey_length)
shapiro.test(df$prey_width)
shapiro.test(df$prey_weight) #it a miracle - no normal distribution

#log transformation of all biometrical data
df_log <- df
df_log["log_prey_length"] <- log10 (df_log$prey_length)
df_log["log_prey_weight"] <- log10 (df_log$prey_weight)
df_log["log_prey_width"] <- log10 (df_log$prey_width)
df_log["log_ant_length"] <- log10 (df_log$ant_length)
df_log["log_ant_width"] <- log10 (df_log$ant_width)
df_log["log_ant_weight"] <- log10 (df_log$ant_weight)

#Checking Normaldistribuiton of log10 data
shapiro.test(df_log$log_ant_length)
shapiro.test(df_log$log_ant_width)
shapiro.test(df_log$log_ant_weight)
shapiro.test(df_log$log_prey_length)
shapiro.test(df_log$log_prey_width)
shapiro.test(df_log$log_prey_weight) #still significant (bc high n) but 0,95 < W < 0,99


# save dataframe as csv for easier exces later on 
write.csv(df_log, "all_ants_log.csv")



#=============================================================================#
#####  Calculating important values and indexes #####
#==============================================================================#

############## Prey Shape (prey length/prey width) 
# index of how long and thin or short and fat prey items are
df["prey_shape"] <- df$prey_length / df$prey_width

############## Prey Area (prey length*prey width) 
# index of how long and thin or short and fat prey items are
df["prey_shape"] <- df$prey_length * df$prey_width

############## Ant Size (ant length*ant width)
# gives the total rectangular 2-dimensional area each ant occupies
df["ant_size"] <- df$ant_length * df$ant_width


#creating new df´s for carrier specific calculatins
df_single <- df %>%
  filter (total_ant_number == 1) #only single carrier

df_multi <- df %>% 
  filter(total_ant_number > 1)  #only multiple carrier
############## Ant loading (prey weight / ant weight)
# Measure of how much an ant carries per unit mass of the ant
# It´s calculated for single and multiple carriers
df_single["ant_loading"] <- df_single$ant_weight / df_single$prey_weight


############## Relative loading single (loading / ant weight)
# relative measure of how heavily loaded ants are relative to their weight.
# single carriers only
df_single["relativ_loading"] <- df_single$ant_loading / (1000*df_single$ant_weight)


############## Load per Ant (loading / ant number)
# mean amount carried per ant for multi-carrier prey items..
# multi carriers only
df_multi["load_per_ant"] <- df_multi$ant_loading / df_multi$total_ant_number

#==============================================================================
                      #### 1. Ant Biometrics Analysis ####
#==============================================================================#

#==============================================================================#
            ##### 1.1 Size differences between Ant species ##### 
#==============================================================================#


#Welch t-test 
t.test(log10(ant_length) ~ ant_species, data = df)
t.test(log10(ant_weight*1000) ~ ant_species, data = df)
t.test(log10(ant_width) ~ ant_species, data = df)

#Effektgröße (Cohen’s d)
cohens_d(log10(ant_length) ~ ant_species, data = df)
cohens_d(log10(ant_weight*1000) ~ ant_species, data = df)
cohens_d(log10(ant_width) ~ ant_species, data = df)


par(mfrow = c(1, 3))
hist(df$ant_length, main = "Ant length", col = "grey",
     xlab = "length [mm]", ylim = c(0, 1000))
hist(df$ant_width,  main = "Ant width [mm]",
     col = "grey",
     xlab = "width [mm]", ylim = c(0, 1000))
hist((1000*df$ant_weight), main = "Ant weight", col = "grey",
     xlab = "Weight [mg]", ylim = c(0, 2500))

#Box plot by species
boxplot(ant_length ~ ant_species, data = df, main = "Ant length by species",
        ylab = "Ant length [mm]", xlab = "")
boxplot(ant_width ~ ant_species, data = df, main = "Ant width by species",
        ylab = "Ant width [mm]", xlab = "")
boxplot((1000*ant_weight)~ ant_species, data = df,
        main = "Ant weight by species", ylab = "Ant weight [mg]", xlab = "")


#==============================================================================#
            ##### 1.2 Allometry between Ants ##### 
#==============================================================================#
###### 1.2.1 Ant length vs. ants width #####

                      # A. Linear Regression Models
#-----------------------------------------------------------------------------#
#base model no main effects
ant_biomet_width1 <- lm(log10(ant_width) ~ log10(ant_length), data = df)

summary(ant_biomet_width1)
confint(ant_biomet_width1)

#model with only main effects
ant_biomet_width2 <- lm(log10(ant_width) ~ log10(ant_length) + ant_species
                       + forest_type , data = df)

summary(ant_biomet_width2)
confint(ant_biomet_width2)

#model with interacting effect specis an main effects 
ant_biomet_width3 <- lm(log10(ant_width) ~ log10(ant_length) * ant_species
                     + forest_type , data = df)

summary(ant_biomet_width3)
confint(ant_biomet_width3)
#plot(rel_load2w)
AIC(ant_biomet_width1, ant_biomet_width2, ant_biomet_width3)
anova(ant_biomet_width1, ant_biomet_width2, ant_biomet_width3)


                      # B. Plot for LRM (ant_biomet_width3)
#-----------------------------------------------------------------------------#


                            # C. Linear Mixed Models
#-----------------------------------------------------------------------------#
#Creating a LMM based on LRM with only main/fixed effects 
#This includes possible random variation due to the raids (adds random effect)
#We do this based on the rusltus of comparison between the models above


ant_biomet_width_lmm <- lmer(log10(ant_width) ~ log10(ant_length) * 
                        ant_species + forest_type  + (1 | raid_ID),
                      data = df)

#checking variation within Raids
VarCorr(ant_biomet_width_lmm)
icc(ant_biomet_width_lmm)
# 17.5% of variance occurs at the raid level → raids are NOT independent

summary(ant_biomet_width_lmm)
AIC(ant_biomet_width1, ant_biomet_width2, 
    ant_biomet_width3,ant_biomet_width_lmm)



                                # D. Plot LMM
#-----------------------------------------------------------------------------#

# Prediction grid
newdat <- expand.grid(
  ant_length    = seq(min(df$ant_length), max(df$ant_length), length = 200),
  ant_species  = levels(df$ant_species),
  forest_type  = "primary",      # choose whichever you want as reference
  raid_ID      = NA              # population-level prediction (random effect = 0)
)

newdat$log_pred <- predict(ant_biomet_width_lmm,
                           newdata = newdat,
                           re.form = NA)

# LMM plot (log-scale axes)
ggplot(df, aes(x = ant_length, y = log10(ant_width), color = ant_species)) +
  geom_point(alpha = 0.35, size = 1) +
  
  # model predicted regression lines
  geom_line(data = newdat,
            aes(x = ant_length, y = log_pred, color = ant_species),
            linewidth = 1.2) +
  
  # log10 axes
  scale_x_log10() +
  scale_y_log10() +
  
  scale_colour_manual(values = c(
    "D. wilverthi" = "#56B4E9",
    "D. sjostedti" = "#E69F00"
  )) +
  
  theme_clean(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 12)
  ) +
  
  labs(
    x = "Ant body length (log10 mm)",
    y = "Ant body width (log10 mm)",
    color = "Species",
    title = "LMM-predicted allometric relationship (log10 scale)"
  )




###### 1.2.2 Ant length vs. ants weight #####
ant_biometrics_2 <- lm(log10(ant_length) ~ log10(ant_weight) * ant_species
                     + forest_type , data = df)

summary(ant_biometrics_2)
confint(ant_biometrics_2)





#=============================================================================
                      #### 2. Prey Biometrics Analysis #### 
#=============================================================================#
##### Preparation of raw data and setting df #####
#reading in data
df_dis <- read.delim(
  file   = "prey_spectra_final.csv",
  sep  = ";",
  dec = "," # <-- very important for comma decimals
)

# Normalize formatting (lowercase + trim white space)
df_dis$dismembered <- tolower(trimws(df_dis$dismembered))

# Keep only valid states ("yes" / "no")
df_dis <- df_dis %>%
  filter(dismembered %in% c("yes", "no"))

table(df_dis$dismembered)

# Make sure its a factor
df_dis$dismembered <- factor(df_dis$dismembered, levels = c("no", "yes"))
df_dis$ant_species <- factor(df_dis$ant_species)
df_dis$forest_type <- factor(df_dis$forest_type)

# Add important factors 
df_dis["prey_shape"] <- df_dis$prey_length / df_dis$prey_width
df_dis["prey_area"] <- df_dis$prey_length * df_dis$prey_width

#==============================================================================
##### 2.1 Size differences in Prey Items ##### 
#==============================================================================#

###### 2.1.1 Testing against Ants #####
#Welch t-test 
t.test(log10(prey_length) ~ ant_species, data = df_dis)
t.test(log10(prey_width) ~ ant_species, data = df_dis)
t.test(log10(prey_weight*1000) ~ ant_species, data = df_dis)
t.test(log10(prey_shape) ~ ant_species, data = df_dis)
t.test(log10(prey_area) ~ ant_species, data = df_dis)

#Effektgröße (Cohen’s d)
cohens_d(log10(prey_length) ~ ant_species, data = df_dis)
cohens_d(log10(prey_width) ~ ant_species, data = df_dis)
cohens_d(log10(prey_weight*1000) ~ ant_species, data = df_dis)
cohens_d(log10(prey_shape) ~ ant_species, data = df_dis)
cohens_d(log10(prey_area) ~ ant_species, data = df_dis)

#Plotting the results - HIER GEHT SCHON WIEDER NIX -BS wrong comand!!
par(mfrow = c(1, 5))
hist(df_dis$prey_length, main = "Prey length", col = "grey",
     xlab = "length [mm]", ylim = c(0, 1000))
hist(df_dis$prey_width,  main = "Prey width [mm]",  col = "grey",
     xlab = "width [mm]", ylim = c(0, 1000))
hist((1000*df_dis$prey_weight), main = "Prey weight", col = "grey",
     xlab = "Weight [mg]", ylim = c(0, 2500))
hist(df_dis$prey_shape, main = "Prey Shape", col = "grey",
     xlab = "Shape (length/width) [mm]", ylim = c(0, 1000))
hist(df_dis$prey_area,  main = "Prey Area [mm]",  col = "grey",
     xlab = "Prey Area [mm^2]", ylim = c(0, 1000))


#Box plot by species
par(mfrow = c(1, 5))
boxplot(prey_length ~ ant_species, data = df_dis,
        main = "Prey length by ant species", ylab = "Prey length [mm]", xlab = "")
boxplot(prey_width ~ ant_species, data = df_dis, 
        main = "Prey width by ant species", ylab = "Prey width [mm]", xlab = "")
boxplot((1000*prey_weight)~ ant_species, data = df_dis,
        main = "Prey weight by species", ylab = "Prey weight [mg]", xlab = "")
boxplot(prey_shape ~ ant_species, data = df_dis, 
        main = "Prey shape by ant species", ylab = "Prey shape", xlab = "")
boxplot(prey_area ~ ant_species, data = df_dis, 
        main = "Prey area by ant species", ylab = "Prey area [mm^2]", xlab = "")


###### 2.1.2 Testing against Forest #####
#Welch t-test 
t.test(log10(prey_length) ~ forest_type, data = df_dis)
t.test(log10(prey_width) ~ forest_type, data = df_dis)
t.test(log10(prey_weight*1000) ~ forest_type, data = df_dis)
t.test(log10(prey_shape) ~ forest_type, data = df_dis)
t.test(log10(prey_area) ~ forest_type, data = df_dis)

#Effektgröße (Cohen’s d)
cohens_d(log10(prey_length) ~ forest_type, data = df_dis)
cohens_d(log10(prey_width) ~ forest_type, data = df_dis)
cohens_d(log10(prey_weight*1000) ~ forest_type, data = df_dis)
cohens_d(log10(prey_shape) ~ forest_type, data = df_dis)
cohens_d(log10(prey_area) ~ forest_type, data = df_dis)


#Plotting the results 
#Box plot by forest
par(mfrow = c(1, 5))
boxplot(prey_length ~ forest_type, data = df_dis,
        main = "Prey length by forest type", ylab = "Prey length [mm]", xlab = "")
boxplot(prey_width ~ forest_type, data = df_dis, 
        main = "Prey width by forest types", ylab = "Prey width [mm]", xlab = "")
boxplot((1000*prey_weight)~ forest_type, data = df_dis,
        main = "Prey weight by forest types", ylab = "Prey weight [mg]", xlab = "")
boxplot(prey_shape ~ forest_type, data = df_dis, 
        main = "Prey shape by forest types", ylab = "Prey shape", xlab = "")
boxplot(prey_area ~ forest_type, data = df_dis, 
        main = "Prey area by forest types", ylab = "Prey area [mm^2]", xlab = "")




#=============================================================================
                      #### 3 Ant-Prey-Relationship  #### 
#=============================================================================#

#=============================================================================#
                     ##### 3.1 Weight-Relationship  #### 
                #Ant weight vs prey weight (single workers) 
#=============================================================================#
###### 3.1.1 Data preparation #####
#setting species and forest as factor (aka categorical variable)
df_single$ant_species <- factor(df_single$ant_species)
df_single$forest_type <- factor(df_single$forest_type)

#In prey_shape and ant_size are INF (when not log transformed data == 1.00)
#We assume the index is just 0 here to get rid of errors 
df_single <- df_single %>%
  mutate(across(c(ant_size, prey_shape),
                ~ ifelse(is.infinite(.), 0, .))) #if Inf/-Inf write 0 

###### 3.1.2 Choosing the fittest Model #####
#base model without explaining variables 
#weight relations for single carriers Nr.1
weight_sc1 <- lm(log10(prey_weight) ~ log10(ant_weight), data = df_single) 

summary(weight_sc1 )
confint(weight_sc1 ) # 95%-Convidence intervall

#model with only main effects
weight_sc2 <- lm(log10(prey_weight) ~ log10(ant_weight) + ant_species + forest_type + 
          prey_shape, data = df_single)

summary(weight_sc2)
confint(weight_sc2)
AIC(weight_sc1 , weight_sc2)  #compares models - Akaike Information Criterion

plot(weight_sc2) #standard plots to check homoskedasticity, linearity, leverage
bptest(weight_sc2) #Breusch-Pagan-test on homoskedasticity (is Varianz in Residuals constant)

#model with interacting effect specis an main effects 
#test this to see if slope is different when species are considered
weight_sc3 <- lm(log10(prey_weight) ~ log10(ant_weight) * ant_species + forest_type + 
           prey_shape, data = df_single)

summary(weight_sc3)
confint(weight_sc3)
AIC(weight_sc1, weight_sc2, weight_sc3)
anova(weight_sc2, weight_sc3)


###### 3.1.3 Plotting: Ant weight vs prey weight (single workers)#####
#This is the plot showing the weight_sc3
#Regression line should be for both species since species has no sig. effect

ggplot(df_single, aes(x = log10(1000*ant_weight), y = log10(1000*prey_weight),
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

#TO DO: Check for used model (is it the right one?)





#=============================================================================#
                    ##### 3.2 Dimensional matching  #### 
                     #Prey size / shape vs Ant size
#=============================================================================#


#=============================================================================#
                   ##### 3.3 Loading in Single Workers  #### 
#=============================================================================#
###### 3.3.1 Relative load vs. Ant weight (single workers) ######
#Relaive load is defined by ant_load / ant_weight

                         # A. Linear Regression Models
#-----------------------------------------------------------------------------#
#base model no main effects
rel_load1w <- lm(log10(relativ_loading) ~ log10(ant_weight), data = df_single) 

summary(rel_load1w)
confint(rel_load1w) # 95%-Convidence intervall
#plot(rel_load1w) 

#model with only main effects
rel_load2w <- lm(log10(relativ_loading) ~ log10(ant_weight) + ant_species +
                   forest_type +  prey_shape, data = df_single)

summary(rel_load2w)
confint(rel_load2w)
#plot(rel_load2w)

#model with interacting effect specis an main effects 
rel_load3w <- lm(log10(relativ_loading) ~ log10(ant_weight) * ant_species + 
                   forest_type + prey_shape, data = df_single)

summary(rel_load3w)
confint(rel_load3w)
AIC(rel_load1w, rel_load2w, rel_load3w)
anova(rel_load1w, rel_load2w, rel_load3w)
#plot(rel_load3w)

                          # B. Linear Mixed Models
#-----------------------------------------------------------------------------#
#Creating a LMM based on LRM with only main/fixed effects 
#This includes possible random variation due to the raids (adds random effect)
#We do this based on the rusltus of comparison between the models above


#The random intercept for raid identity explained 13% of the variance (ICC=0.13)

rel_load_LMMw <- lmer(log10(relativ_loading) ~ log10(ant_weight) + 
                       ant_species + forest_type + prey_shape + (1 | raid_ID),
                 data = df_single)

summary(rel_load_LMMw)
AIC(rel_load1w, rel_load2w, rel_load3w, rel_load_LMMw)

#checking variation within Raids
VarCorr(rel_load_LMMw)
icc(rel_load_LMMw)
                                  # C. DHARMa 
#-----------------------------------------------------------------------------#

sim_rel_load_LMMw <- simulateResiduals(rel_load_LMMw)  #Simulating Residuals

# plot(sim_rel_load_LMMw ) #plot diagnostic

# Testing Distribution & Dispersion
testDispersion(sim_rel_load_LMMw )
testUniformity(sim_rel_load_LMMw )
testResiduals(sim_rel_load_LMMw )  



                              # C. Plot for LMM 
#-----------------------------------------------------------------------------#
# Sequencing Ant weight to have fixed points for the predict grid
ant_seq_rel_loadw <- seq(
  from = min(df_single$ant_weight, na.rm = TRUE),
  to   = max(df_single$ant_weight, na.rm = TRUE),
  length.out = 100
)

#Prediction-Grid (for both spices, primary forest, mean prey_shape)
#To calculate predicts and generate one line all effects net to be standardized
newdat_rel_loadw <- expand.grid(
  ant_weight   = ant_seq_rel_loadw,
  ant_species  = levels(df_single$ant_species),
  forest_type  = "primary",
  prey_shape   = mean(df_single$prey_shape, na.rm = TRUE)
)

#LMM-Predicts (only Fixed Effects)

newdat_rel_loadw$pred_log_rel_load <- predict(
  rel_load_LMMw,
  newdata = newdat_rel_loadw,
  re.form = NA   # no random effects in line
)

#Generating Plot with data points + LMM-Lines
plot_rel_weight<- ggplot(df_single, aes(x = log10(ant_weight),
                      y = log10(relativ_loading),
                      colour = ant_species)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_line(data = newdat_rel_loadw,
            aes(x = log10(ant_weight),
                y = pred_log_rel_load,
                colour = ant_species),
            linewidth = 1.1) +
  theme_classic(base_size = 13) +
  labs(
    x = "Ant weight (log10 mg)",
    y = "Relative load (log10 prey mass / ant mass)",
    colour = "Ant species",
    title = "Relative load vs. ant weight (single carriers, LMM predictions)"
  )

# file save
jpeg(file = "Relative Load vs. Ant Weight.jpg",
     width = 25, height = 18, units = "cm", res = 300)
plot_rel_weight
dev.off()


#_____________________________________________________________________________#

###### 3.3.2 Relative load vs. Ant size (single workers) ######

# Relative load is defined by ant_load / ant-weight

                       # A. Linear Regression Models
#-----------------------------------------------------------------------------#
# base model, no main effects
rel_load1s <- lm(log10(relativ_loading) ~ log10(ant_size), data = df_single) 

summary(rel_load1s)
confint(rel_load1s)

# model with only main effects
rel_load2s <- lm(log10(relativ_loading) ~ log10(ant_size) + ant_species +
                   forest_type + prey_shape, data = df_single)

summary(rel_load2s)
confint(rel_load2s)

# model with interaction species × ant_size
rel_load3s <- lm(log10(relativ_loading) ~ log10(ant_size) * ant_species + 
                   forest_type + prey_shape, data = df_single)

summary(rel_load3s)
confint(rel_load3s)
AIC(rel_load1s, rel_load2s, rel_load3s)
anova(rel_load1s, rel_load2s, rel_load3s)


                          # B. Linear Mixed Model
#-----------------------------------------------------------------------------#
# LMM based on best LRM structure, including raid-level random intercept

rel_load_LMMs <- lmer(log10(relativ_loading) ~ log10(ant_size) + 
                        ant_species + forest_type + prey_shape + (1 | raid_ID),
                      data = df_single)

summary(rel_load_LMMs)
AIC(rel_load1s, rel_load2s, rel_load3s, rel_load_LMMs)

# random-effect variance
VarCorr(rel_load_LMMs)
icc(rel_load_LMMs)


                                # C. DHARMa 
#-----------------------------------------------------------------------------#

sim_rel_load_LMMs <- simulateResiduals(rel_load_LMMs)

# plot(sim_rel_load_LMMs)

testDispersion(sim_rel_load_LMMs)
testUniformity(sim_rel_load_LMMs)
testResiduals(sim_rel_load_LMMs)


                              # D. Plot for LMM 
#-----------------------------------------------------------------------------#
# 1) Create ant_size sequence for prediction grid
ant_seq_rel_loads <- seq(
  from = min(df_single$ant_size, na.rm = TRUE),
  to   = max(df_single$ant_size, na.rm = TRUE),
  length.out = 100
)

# 2) Prediction grid (both species, primary forest, mean prey shape)
newdat_rel_loads <- expand.grid(
  ant_size     = ant_seq_rel_loads,
  ant_species  = levels(df_single$ant_species),
  forest_type  = "primary",
  prey_shape   = mean(df_single$prey_shape, na.rm = TRUE)
)

# 3) Fixed-effect LMM predictions
newdat_rel_loads$pred_log_rel_load <- predict(
  rel_load_LMMs,
  newdata = newdat_rel_loads,
  re.form = NA
)

# 4) Create plot
plot_rel_size <- ggplot(df_single, aes(x = log10(ant_size),
                                       y = log10(relativ_loading),
                                       colour = ant_species)) +
  geom_point(alpha = 0.3, size = 1) +
  geom_line(data = newdat_rel_loads,
            aes(x = log10(ant_size),
                y = pred_log_rel_load,
                colour = ant_species),
            linewidth = 1.1) +
  theme_classic(base_size = 13) +
  labs(
    x = "Ant size (log10 mm²)",
    y = "Relative load (log10 prey mass / ant mass)",
    colour = "Ant species",
    title = "Relative load vs. ant size (single carriers, LMM predictions)"
  )

# Save file
jpeg(file = "Relative Load vs. Ant Size.jpg",
     width = 25, height = 18, units = "cm", res = 300)
plot_rel_size
dev.off()

###############################################################################
#Load per ant vs ant weight + single vs multiple carriers######################
#setting species and forest as factor (aka categorical variable)
df_multi$ant_species <- factor(df_multi$ant_species)
df_multi$forest_type <- factor(df_multi$forest_type)

#In prey_shape and ant_size are INF (when not log transformed data == 1.00)
#We assume the index is just 0 here to get rid of errors 
df_multi <- df_multi %>%
  mutate(across(c(ant_size, prey_shape),
                ~ ifelse(is.infinite(.), 0, .))) #if Inf/-Inf write 0 




