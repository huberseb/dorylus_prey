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
df["prey_shape"] <- df$log_prey_length / df$log_prey_width

############## Ant Size (ant length*ant width)
# gives the total rectangular 2-dimensional area each ant occupies
df["ant_size"] <- df$log_ant_length / df$log_ant_width

############## Ant loading (prey weight / ant weight)
# Measure of how much an ant carries per unit mass of the ant
# It´s calculated for single and multiple carriers
df["ant_loading"] <- df$log_prey_weight / df$log_ant_weight


#creating new df´s for carrier specific calculatins
df_single <- df %>%
  filter (total_ant_number == 1) #only single carrier

df_multi <- df %>% 
  filter(total_ant_number > 1)  #only multiple carrier


############## Relative loading single (loading / ant weight)
# relative measure of how heavily loaded ants are relative to their weight.
# single carriers only
df_single["relativ_loading"] <- df_single$ant_loading / df_single$log_ant_weight


############## Relative loading multi (loading / ant weight)
# relative measure of how heavily loaded ants are relative to their weight.
# multi carriers only
df_multi["realative_loading"] <- df_multi$ant_loading / df_multi$log_ant_weight



#==============================================================================
                #### Size differences between ant species #### 
#==============================================================================


#Welch t-test 
t.test(ant_length ~ ant_species, data = df)

#Effektgröße (Cohen’s d)

cohens_d(ant_length ~ ant_species, data = df)

par(mfrow = c(1, 3))
hist(df$ant_length, main = "Length", col = "grey")
hist(df$ant_width,  main = "Width",  col = "grey")
hist(df$ant_weight, main = "Weight", col = "grey")

#Box plot by species
boxplot(ant_length ~ ant_species, data = df, main = "Ant length by species")
boxplot(ant_width ~ ant_species, data = df, main = "Ant width by species")
boxplot(ant_weight~ ant_species, data = df, main = "Ant weight by species")



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


#Ant weight vs prey weight (single workers)###################################
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


