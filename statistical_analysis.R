library(ggplot2)
library(effectsize)

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

####df might not set right. TO DO SET DF AND STRUCTURE FROM HERRE ON 

              #Prey Shape (prey length/prey width) 
# index of how long and thin or short and fat prey items are

df_log["prey_shape"] <- df$prey_length / df$prey_width

              #Ant Size (ant length*ant width)
# gives the total rectangular 2-dimensional area each ant occupies
df_index["ant_size"] <- df$ant_length / df$ant_width

              #Ant loading (prey weight / ant weight)
# Measure of how much an ant carries per unit mass of the ant
# It´s calculated for single and multiple carriers
df_index["ant_loading"] <- df$prey_weight / df$ant_weight

              #Relative loading



#### Size differences between ant species ####

# Welch t-test (Standard in R)
t.test(ant_length ~ ant_species, data = df)

# Effektgröße (Cohen’s d)

cohens_d(ant_length ~ ant_species, data = df)

par(mfrow = c(1, 3))
hist(df$ant_length, main = "Length", col = "grey")
hist(df$ant_width,  main = "Width",  col = "grey")
hist(df$ant_weight, main = "Weight", col = "grey")

# oder Boxplots nach Art
boxplot(ant_length ~ ant_species, data = df, main = "Ant length by species")
boxplot(ant_width ~ ant_species, data = df, main = "Ant width by species")
boxplot(ant_weight~ ant_species, data = df, main = "Ant weight by species")