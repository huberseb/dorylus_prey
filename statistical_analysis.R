library(ggplot2)
library(effectsize)
library(tidyverse)
library(dplyr)
library(lmtest)
library(performance)

#reading in data
df <- read.delim(
  file   = "all_ants_all_prey.csv",
  sep  = ",",
  dec = "." # <-- very important for comma decimals
)

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
#reading in data
df <- read.delim(
  file   = "all_ants_log.csv",
  sep  = ",",
  dec = "." # <-- very important for comma decimals
) 


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
hist(df$ant_width,  main = "Ant width [mm]",  col = "grey",
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
ant_biometrics_1 <- lm(log10(ant_length) ~ log10(ant_width) * ant_species
                     + forest_type , data = df)

summary(ant_biometrics_1)
confint(ant_biometrics_1)

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
                     #### 3.1 Weight-Relationship  #### 
#=============================================================================#

#setting species and forest as factor (aka categorical variable)
df_single$ant_species <- factor(df_single$ant_species)
df_single$forest_type <- factor(df_single$forest_type)

#In prey_shape and ant_size are INF (when not log transformed data == 1.00)
#We assume the index is just 0 here to get rid of errors 
df_single <- df_single %>%
  mutate(across(c(ant_size, prey_shape),
                ~ ifelse(is.infinite(.), 0, .))) #if Inf/-Inf write 0 


#Ant weight vs prey weight (single workers)#####################################
#base model without explaining variables 
#weight relations for single carriers Nr.1
weight_sc1 <- lm(log_prey_weight ~ log_ant_weight, data = df_single) 

summary(weight_sc1 )
confint(weight_sc1 ) # 95%-Convidence intervall

#model with only main effects
weight_sc2 <- lm(log_prey_weight ~ log_ant_weight + ant_species + forest_type + 
          prey_shape, data = df_single)

summary(weight_sc2)
confint(weight_sc2)
AIC(weight_sc1 , weight_sc2)  #compares models - Akaike Information Criterion

plot(weight_sc2) #standard plots to check homoskedasticity, linearity, leverage
bptest(weight_sc2) #Breusch-Pagan-test on homoskedasticity (is Varianz in Residuals constant)

#model with interacting effect specis an main effects 
#test this to see if slope is different when species are considered
weight_sc3 <- lm(log_prey_weight ~ log_ant_weight * ant_species + forest_type + 
           prey_shape, data = df_single)

summary(weight_sc3)
confint(weight_sc3)
AIC(weight_sc1, weight_sc2, weight_sc3)
anova(weight_sc2, weight_sc3)

#Relative load vs. Ant weight (single workers)################################
#base model no main effects
rel_load1 <- lm(log10(relativ_loading) ~ log10(ant_weight), data = df_single) 

summary(rel_load1)
confint(rel_load1) # 95%-Convidence intervall
#plot(rel_load1) 

#model with only main effects
rel_load2 <- lm(log10(relativ_loading) ~ log10(ant_size) + ant_species + forest_type + 
            prey_shape, data = df_single)

summary(rel_load2)
confint(rel_load2)
#plot(rel_load2)

#model with interacting effect specis an main effects 
rel_load3 <- lm(log10(relativ_loading) ~ log10(ant_size) * ant_species + forest_type + 
            prey_shape, data = df_single)

summary(rel_load3)
confint(rel_load3)
AIC(rel_load1, rel_load2, rel_load3)
anova(rel_load1, rel_load2, rel_load3)
#plot(rel_load3)


#Load per ant vs ant weight + single vs multiple carriers######################
#setting species and forest as factor (aka categorical variable)
df_multi$ant_species <- factor(df_multi$ant_species)
df_multi$forest_type <- factor(df_multi$forest_type)

#In prey_shape and ant_size are INF (when not log transformed data == 1.00)
#We assume the index is just 0 here to get rid of errors 
df_multi <- df_multi %>%
  mutate(across(c(ant_size, prey_shape),
                ~ ifelse(is.infinite(.), 0, .))) #if Inf/-Inf write 0 




