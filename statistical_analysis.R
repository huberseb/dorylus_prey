library(ggplot2)
library(effectsize)
library(tidyvers)
library(dplyr)
library(lmtest)
library(performance)

#reading in data
df <- read.delim(
  file   = "all_ants_all_prey.csv",
  sep  = ",",
  dec = "." # <-- very important for comma decimals
)

#==============================================================================
                      #### Testing Distribution ####
#==============================================================================

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



#==============================================================================
            #### Calculating important values and indexes ####
#==============================================================================
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

############## Ant loading (prey weight / ant weight)
# Measure of how much an ant carries per unit mass of the ant
# It´s calculated for single and multiple carriers
df_single["ant_loading"] <- df_single$ant_weight / df_single$prey_weight


#creating new df´s for carrier specific calculatins
df_single <- df %>%
  filter (total_ant_number == 1) #only single carrier

df_multi <- df %>% 
  filter(total_ant_number > 1)  #only multiple carrier


############## Relative loading single (loading / ant weight)
# relative measure of how heavily loaded ants are relative to their weight.
# single carriers only
df_single["relativ_loading"] <- df_single$ant_loading / df_single$ant_weight


############## Load per Ant (loading / ant number)
# mean amount carried per ant for multi-carrier prey items..
# multi carriers only
df_multi["load_per_ant"] <- df_multi$ant_loading / df_multi$total_ant_number


#==============================================================================#
                #### Size differences between ant species #### 
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
hist(df$ant_length, main = "Ant length", col = "grey", xlab = "length [mm]", ylim = c(0, 1000))
hist(df$ant_width,  main = "Ant width [mm]",  col = "grey", xlab = "width [mm]", ylim = c(0, 1000))
hist((1000*df$ant_weight), main = "Ant weight", col = "grey", xlab = "Weight [mg]", ylim = c(0, 2500))

#Box plot by species
boxplot(ant_length ~ ant_species, data = df, main = "Ant length by species",
        ylab = "Ant length [mm]", xlab = "")
boxplot(ant_width ~ ant_species, data = df, main = "Ant width by species",
        ylab = "Ant width [mm]", xlab = "")
boxplot((1000*ant_weight)~ ant_species, data = df, main = "Ant weight by species",
        ylab = "Ant weight [mg]", xlab = "")




#==============================================================================
                     #### Linear Regression models #### 
#==============================================================================

#setting species and forest as factor (aka categorical variable)#
df_single$ant_species <- factor(df_single$ant_species)
df_single$forest_type <- factor(df_single$forest_type)

#In prey_shape and ant_size are INF (when not log transformed data == 1.00)
#We assume the index is just 0 here to get rid of errors 
df_single <- df_single %>%
  mutate(across(c(ant_size, prey_shape),
                ~ ifelse(is.infinite(.), 0, .))) #if Inf/-Inf write 0 


#Ant weight vs prey weight (single workers)#####################################
#base model without explaining variables
lm1 <- lm(log_prey_weight ~ log_ant_weight, data = df_single) 

summary(lm1)
confint(lm1) # 95%-Convidence intervall

#model with only main effects
lm2 <- lm(log_prey_weight ~ log_ant_weight + ant_species + forest_type + 
          prey_shape, data = df_single)

summary(lm2)
confint(lm2)
AIC(lm1, lm2)  #compares models - Akaike Information Criterion

plot(lm2) #standard plots to check homoskedasticity, linearity, leverage
bptest(lm2) #Breusch-Pagan-test on homoskedasticity (is Varianz in Residuals constant)

#model with interacting effect specis an main effects 
#test this to see if slope is different when species are considered
lm3 <- lm(log_prey_weight ~ log_ant_weight * ant_species + forest_type + 
           prey_shape, data = df_single)

summary(lm3)
confint(lm3)
AIC(lm1, lm2, lm3)
anova(lm2, lm3)

#Relative load vs. Ant weight (single workers)################################
#base model no main effects
lm4 <- lm(log10(relativ_loading) ~ log10(ant_weight), data = df_single) 

summary(lm4)
confint(lm4) # 95%-Convidence intervall
plot(lm4) 

#model with only main effects
lm5 <- lm(log10(relativ_loading) ~ log10(ant_size) + ant_species + forest_type + 
            prey_shape, data = df_single)

summary(lm5)
confint(lm5)
plot(lm5)

#model with interacting effect specis an main effects 
lm6 <- lm(log10(relativ_loading) ~ log10(ant_size) * ant_species + forest_type + 
            prey_shape, data = df_single)

summary(lm6)
confint(lm6)
AIC(lm4, lm5, lm6)
anova(lm4, lm5, lm6)
plot(lm6)


#Load per ant vs ant weight + single vs multiple carriers######################
#setting species and forest as factor (aka categorical variable)
df_multi$ant_species <- factor(df_multi$ant_species)
df_multi$forest_type <- factor(df_multi$forest_type)

#In prey_shape and ant_size are INF (when not log transformed data == 1.00)
#We assume the index is just 0 here to get rid of errors 
df_multi <- df_multi %>%
  mutate(across(c(ant_size, prey_shape),
                ~ ifelse(is.infinite(.), 0, .))) #if Inf/-Inf write 0 




