library(ggplot2)                                  #        ,_     _,
library(ggthemes)                                 #          '._.'
library(patchwork)                                #     '-,   (_)   ,-'     D
library(effectsize)                               #       '._ .:. _.'       O
library(tidyverse)                                #        _ '|Y|' _        R 
library(dplyr)                                    #      ,` `>\ /<` `,      Y
library(lmtest)                                   #     ` ,-`  I  `-, `     L
library(performance)  #e.g. for ICC calculations  #       |   /=\   |       U
library(lme4)         #for LMM´s                  #     ,-'   |=|   '-,     S
library(DHARMa)                                   #          ( = )        
library(ggpattern)                                #           \_/       
library(patchwork)   #ombining plots              # # # # # # # # # # # # # # #  
library(cowplot) 
library(gghalves)    #violin plots in boxplots
library(ggforce)
library(kSamples)    #KS and DA Test (cite)
library(diptest)      #Dip test (cite)
library(multimode)   #Silverman´s test (cite)
library(ggbeeswarm)  #density bases point clouds -not necce
library(car)         #Varianz test 
library(flextable)   #Tabel tool
library(officer)     # creates .docx files
library(lmtest)      # BP test (cite)
library(sandwich)    # HC3 p values (cite)

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
df["prey_area"] <- df$prey_length * df$prey_width

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

#===================================.==========================================

                    #### 1. Ant Biometrics Analysis ####

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

#Kolmogorov–Smirnov-Test
#compares whole distribution between two groups
#non parametric and sensitive vor varianzes 
#empfindlich auf globalen Unterschied

#Pull biometric data per species
len_wil <- df %>% filter(ant_species == "D. wilverthi") %>% pull(ant_length)
len_sjo <- df %>% filter(ant_species == "D. sjostedti") %>% pull(ant_length)
wid_wil <- df %>% filter(ant_species == "D. wilverthi") %>% pull(ant_width)
wid_sjo <- df %>% filter(ant_species == "D. sjostedti") %>% pull(ant_width)
wei_wil <-df%>% filter(ant_species == "D. wilverthi") %>% pull(ant_weight)*1000
wei_sjo <-df%>% filter(ant_species == "D. sjostedti") %>% pull(ant_weight)*1000


#K-S Test
ks_len <- ks.test(len_wil, len_sjo)
ks_len
ks_wid <- ks.test(wid_wil, wid_sjo)
ks_wid
ks_wei <- ks.test(wei_wil, wei_sjo)
ks_wei


# Anderson–Darling-Test
#more sensitive in the taisl, then KS
#empfindlich auf Tail-Unterschiede
ad_len <- ad.test(len_wil, len_sjo)
ad_len
ad_wid <- ad.test(wid_wil, wid_sjo)
ad_wid
ad_wei <- ad.test(wei_wil, wei_sjo)
ad_wei

##------------------ Hartigan´s dip test for multimodality ------------------##
# here we test if there are more than one peak in our data 

# Length: 
dip_len_wil <- dip.test(len_wil)   
dip_len_sjo <- dip.test(len_sjo)  

# Width: 
dip_wid_wil <- dip.test(wid_wil)  
dip_wid_sjo <- dip.test(wid_sjo)  

# Weight
dip_wei_wil <- dip.test(wei_wil)  
dip_wei_sjo <- dip.test(wei_sjo)

# Rusults
dip_len_wil; dip_len_sjo
dip_wid_wil; dip_wid_sjo
dip_wei_wil; dip_wei_sjo



##------------------ Silverman's test for number of modes ------------------##
#here we test if the distribution has a maximum of k modi (peaks)
#Set R for bootstraping to 999 for easy computing. 

# D. wilverthi – length
silv_len_wil_k1 <- modetest(len_wil, mod0 = 1, method = "SI", B = 999)  
silv_len_wil_k2 <- modetest(len_wil, mod0 = 2, method = "SI", B = 999)
silv_len_wil_k3 <- modetest(len_wil, mod0 = 3, method = "SI", B = 999)

# D. sjostedti – length
silv_len_sjo_k1 <- modetest(len_sjo, mod0 = 1, method = "SI", B = 999)
silv_len_sjo_k2 <- modetest(len_sjo, mod0 = 2, method = "SI", B = 999)
silv_len_sjo_k3 <- modetest(len_sjo, mod0 = 3, method = "SI", B = 999)

# D. wilverthi – width
silv_wid_wil_k1 <- modetest(wid_wil, mod0 = 1, method = "SI", B = 999)
silv_wid_wil_k2 <- modetest(wid_wil, mod0 = 2, method = "SI", B = 999)
silv_wid_wil_k3 <- modetest(wid_wil, mod0 = 3, method = "SI", B = 999)

# D. sjostedti – width
silv_wid_sjo_k1 <- modetest(wid_sjo, mod0 = 1, method = "SI", B = 999)
silv_wid_sjo_k2 <- modetest(wid_sjo, mod0 = 2, method = "SI", B = 999)
silv_wid_sjo_k3 <- modetest(wid_sjo, mod0 = 3, method = "SI", B = 999)

# D. wilverthi – weight (mg)
silv_wei_wil_k1 <- modetest(wei_wil, mod0 = 1, method = "SI", B = 999)
silv_wei_wil_k2 <- modetest(wei_wil, mod0 = 2, method = "SI", B = 999)
silv_wei_wil_k3 <- modetest(wei_wil, mod0 = 3, method = "SI", B = 999)

# D. sjostedti – weight (mg)
silv_wei_sjo_k1 <- modetest(wei_sjo, mod0 = 1, method = "SI", B = 999)
silv_wei_sjo_k2 <- modetest(wei_sjo, mod0 = 2, method = "SI", B = 999)
silv_wei_sjo_k3 <- modetest(wei_sjo, mod0 = 3, method = "SI", B = 999)

# Results
silv_len_wil_k1; silv_len_wil_k2; silv_len_wil_k3
silv_len_sjo_k1; silv_len_sjo_k2; silv_len_sjo_k3
silv_wid_wil_k1; silv_wid_wil_k2; silv_wid_wil_k3
silv_wid_sjo_k1; silv_wid_sjo_k2; silv_wid_sjo_k3
silv_wei_wil_k1; silv_wei_wil_k2; silv_wei_wil_k3
silv_wei_sjo_k1; silv_wei_sjo_k2; silv_wei_sjo_k3




            ######### A. Histograms #######
#-----------------------------------------------------------------------------#
#overview
par(mfrow = c(1, 3))
hist(df$ant_length, main = "Ant length", col = "grey",
     xlab = "length [mm]", ylim = c(0, 1000))
hist(df$ant_width,  main = "Ant width [mm]",
     col = "grey",
     xlab = "width [mm]", ylim = c(0, 1000))
hist((1000*df$ant_weight), main = "Ant weight", col = "grey",
     xlab = "Weight [mg]", ylim = c(0, 2500))

#Detailed Histograms for each metric showing Proportions

#gettin the n for ant species
N_wil <- sum(df$ant_species == "D. wilverthi", na.rm = TRUE)
N_sjo <- sum(df$ant_species == "D. sjostedti", na.rm = TRUE)

#In geom_histo we then give every sample a weight: 1/N_species
#..count.. summs up all weights in bin -> gives proportion per species 

#Defining colour
ant_colour <- c ("D. wilverthi" = "#56B4E9", "D. sjostedti" = "#E69F00")

dorylus_colour <- list (
  scale_fill_manual(values = ant_colour),
  scale_colour_manual(values = ant_colour))

#Defining plot style
plot_style <- list(
   dorylus_colour,
  guides(colour = "none"),
  theme_clean(base_size = 18),
  theme(plot.title = element_text(hjust = 1)
  ))


                             #-----LENGTH-----#
plot_ant_length <- ggplot(df, aes(x = ant_length,
                             colour = ant_species, fill =  ant_species)) +
  geom_histogram( 
    aes( weight = ifelse(ant_species == "D. wilverthi", 1 / N_wil, 1 / N_sjo),
      y = after_stat(..count..) ),  
    position = "identity",
    alpha    = 0.5,
    bins     = 30) +
  plot_style + 
  labs(
    x = "Ant body length [mm]",
    y = "Proportion per species",
    fill = "Dorylus species",
  )
plot_ant_length

                             #-----WIDTH-----#
plot_ant_width <- ggplot(df, aes(x = ant_width,
                             colour = ant_species, fill =  ant_species)) +
  geom_histogram( 
    aes( weight = ifelse(ant_species == "D. wilverthi", 1 / N_wil, 1 / N_sjo),
         y = after_stat(..count..) ),  
    position = "identity",
    alpha    = 0.5,
    bins     = 30) +
  plot_style +  
  labs(
    x = "Ant head width [mm]",
    y = "Proportion per species",
    fill = "Dorylus species",
  )
plot_ant_width

                            #-----WEIGHT-----#
plot_ant_weight <- ggplot(df, aes(x = (1000*ant_weight),
                                 colour = ant_species, fill =  ant_species)) +
  geom_histogram( 
    aes( weight = ifelse(ant_species == "D. wilverthi", 1 / N_wil, 1 / N_sjo),
         y = after_stat(..count..) ),  
    position = "identity",
    alpha    = 0.5,
    bins     = 30) +
  plot_style +  
  labs(
    x = "Ant weight [mg]",
    y = "Proportion per species",
    fill = "Dorylus species",
  )
plot_ant_weight



            ######### B. Combining all Histogramms#######

#defining universal legend theme 
ant_legend_theme <- theme(
  legend.position      = c(0.98, 0.99),           # inwards, top right
  legend.justification = c("right", "top"),       # anker point
  legend.background    = element_rect(fill = "white", colour = NA),
  legend.box.margin    = margin(0, 0, 0, 0),
  legend.margin        = margin(0, 0, 0, 0),
  legend.key.size      = unit(0.7, "lines"),      # compact legend
  legend.key.spacing   = unit(0.5, "lines"),
  legend.key.spacing.x = unit(1.2, "lines"),
  legend.text          = element_text(size = 18, face = "italic"),
  legend.title         = NULL, 
  legend.text.align  = 0.5,
  legend.title.align = 0.5
)

#defining clean panel theme 
clean_panel_theme <- theme(
  axis.title.x = element_text(size = 18, margin = margin(t = 10)),
  axis.title.y = element_text(size = 18, margin = margin(r = 14)),
  axis.text.x  = element_text(size = 14),
  axis.text.y  = element_text(size = 14),
  panel.border      = element_blank(),     # kein Rahmen um das Panel
  panel.grid.major  = element_blank(),     # keine gestrichelten Linien
  panel.grid.minor  = element_blank(),
  panel.background  = element_blank(),     # keine graue Fläche
  plot.background   = element_blank(),
)

#defining scale range 
#lims_length <- range(df$ant_length, na.rm = TRUE)
#Ant lengthe was a bit more tricky needed workaround
pb <- ggplot_build(plot_ant_length)
xlims   <- pb$layout$panel_params[[1]]$x.range
xbreaks <- pb$layout$panel_params[[1]]$x$breaks

scale_length <- scale_x_continuous( limits = xlims, breaks = xbreaks)
lims_weight <- c( -1, 44)
lims_width <- c(0.5,4.3)

#Adding all themes to plots
#adding y title, overall panel theme and legend 
plot_ant_length  <- plot_ant_length  + 
  labs(y = "Proportion per species") +
  clean_panel_theme + theme(legend.position = "none") + scale_length

plot_ant_width   <- plot_ant_width   + 
  labs(y = NULL) + clean_panel_theme +
  theme(legend.position = "none") +
  scale_x_continuous(limits = lims_width, breaks = pretty(lims_width, n = 4))

plot_ant_weight  <- plot_ant_weight  + labs(y = NULL) + clean_panel_theme +
  ant_legend_theme + 
  scale_x_continuous(limits = lims_weight, breaks = pretty(lims_weight, n = 5))


#Combining everything
plot_ant_biometrics <- (plot_ant_length + plot_ant_width + plot_ant_weight) +
  theme(
    plot.margin      = margin(5, 20, 5, 20),
    panel.spacing    = unit(0.8, "lines")  # Abstand zwischen Panels
  )

# file save
jpeg(file = "Histogramm Biometrics Ants.jpg",
     width = 60, height = 18, units = "cm", res = 300)
plot_ant_biometrics
dev.off()





            ######### C. Boxplots #######
#Overview
par(mfrow = c(1, 3))
boxplot(ant_length ~ ant_species, data = df, main = "Ant length by species",
        ylab = "Ant length [mm]", xlab = "")
boxplot(ant_width ~ ant_species, data = df, main = "Ant width by species",
        ylab = "Ant width [mm]", xlab = "")
boxplot((1000*ant_weight)~ ant_species, data = df,
        main = "Ant weight by species", ylab = "Ant weight [mg]", xlab = "")

#Universal Boxplot esthetics
boxp_style <- list(
  geom_boxplot( 
  width         = 0.5,     # width of boxes 
  size          = 0.6,# size of lines
  colour        = "grey11",
  alpha         = 0.55, 
  outliers      = FALSE, 
  outlier.shape = NA,
  outlier.colour = "grey11",
  outlier.fill   = NA,
  outlier.size  = 1,
  outlier.stroke = 0.3,
  outlier.alpha = 0.3, 
  staplewidth   = 0.5),
#    theme(panel.grid = element_blank(),
 #         axis.ticks.x = element_blank(),  
  #        axis.ticks.length = unit(0, "pt")),
  dorylus_colour,
  guides(colour = "none"),
  theme_clean(base_size = 18) %+replace% theme(
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_blank(),
    panel.border     = element_blank(),
    plot.title = element_text(hjust = 1)          # remove tick length
  ))


#Point cloud with outliners 
boxp_distribution <- list( 
  geom_sina(
    size = 1,
    alpha = 0.3,
    color = "grey82",
    method = "density",
    position = position_nudge(x = 0)
  ))

#Creaiton violin overlay   NOTE: for now disabled in idv. plots
boxp_violin <- list(
  geom_half_violin(
    side        = "r",      # right side
    width       = 0.7,      # width from center
    trim        = FALSE,    # show whole distribution
    alpha       = 0.4,     
    colour      = NA,       # no rim 
    #fill        = an_species,
    position    = position_nudge(0.08, 0) 
  ))

#losing all x axis and grids
without_x_achse <-   theme(
  axis.title.x = element_blank(),
  axis.text.x  = element_blank(),
  axis.ticks.x = element_blank() )

#losing all y axis
without_y_achse <-  theme(
  axis.title.y = element_blank(),
  axis.text.y  = element_blank(),
  axis.ticks.y = element_blank() )

#Legend for Boxplots
boxp_legend_theme <- theme(
  legend.position      = c(0.98, 0.5),           # inwards, top right
  legend.justification = c("right", "center"),       # anker point
  legend.background    = element_rect(fill = "white", colour = NA),
  legend.box.margin    = margin(0, 0, 0, 0),
  legend.margin        = margin(0, 0, 0, 0),
  legend.key.size      = unit(0.7, "lines"),      # compact legend
  legend.key.spacing   = unit(0.5, "lines"),
  legend.key.spacing.x = unit(1.2, "lines"),
  legend.text          = element_text(size = 18, face = "italic"),
  legend.title         =  NULL, 
  legend.text.align  = 0.5,
  legend.title.align = 0.5
)

#Length
boxplot_ant_length <- ggplot(df, aes(x = ant_species, y = (ant_length), 
                                colour = ant_species, fill =  ant_species)) +
  boxp_distribution +
  boxp_style +
  #boxp_violin +
  labs( y = "Ant body length [mm]",
        x ="",
        fill = NULL)+
  without_y_achse +
  #without_x_achse +
  coord_flip() +
  boxp_legend_theme 
boxplot_ant_length


#Width
boxplot_ant_width <- ggplot(df, aes(x = ant_species, y = (ant_width), 
                                     colour = ant_species, fill =  ant_species)) +
  boxp_distribution +
  boxp_style +
  #boxp_violin +
  labs( y = "Ant body width [mm]",
        x ="",
        fill = NULL)+
  without_y_achse +
  #without_x_achse +
  coord_flip() +
  boxp_legend_theme 
boxplot_ant_width

#Weight
boxplot_ant_weight <- ggplot(df, aes(x = ant_species, y = (1000*ant_weight), 
                                  colour = ant_species, fill =  ant_species)) +
  boxp_distribution +
  boxp_style +
  #boxp_violin +
  labs( y = "Ant weight [mg]",
        x ="",
        fill = NULL)+
  without_y_achse +
  #without_x_achse +
  coord_flip() +
  boxp_legend_theme 
boxplot_ant_weight




            ######### D. Combining all Boxplots #######

#Panel definition
boxp_panel <- theme( 
  panel.border      = element_blank(),     # no frame around panel
  panel.grid.major  = element_blank(),     # no doted lines
  panel.grid.minor  = element_blank(),
  panel.background  = element_blank(),     
  plot.background   = element_blank(),
  )

#Preperation of individual plots
boxplot_ant_length_a  <- boxplot_ant_length + boxp_panel +
  theme(legend.position = "none") 
boxplot_ant_width_a   <- boxplot_ant_width + boxp_panel +
  theme(legend.position = "none") 
boxplot_ant_weight_a  <- boxplot_ant_weight + boxp_legend_theme + boxp_panel

#Combining everything
boxplot_ant_biometrics <- (boxplot_ant_length_a + boxplot_ant_width_a + 
                             boxplot_ant_weight_a) +
  theme( plot.margin      = margin(5, 20, 5, 20), # between Panels
         panel.spacing    = unit(0.8, "lines"),
         )      

# file save
jpeg(file = "Boxplot Biometrics Ants_v2.jpg",
     width = 60, height = 18, units = "cm", res = 300)
boxplot_ant_biometrics
dev.off()

            ######### E. Combining all Plots #######

#Preparing Boxplots
boxplot_ant_length_c  <- boxplot_ant_length_a  +  labs (y = "") +
  scale_y_continuous(position = "right", limits = xlims, breaks = xbreaks)
boxplot_ant_width_c   <- boxplot_ant_width_a   + labs (y = "") +
  scale_y_continuous(position = "right", limits = lims_width,
                     breaks = pretty(lims_width, n = 4))
boxplot_ant_weight_c  <- boxplot_ant_weight_a + theme(legend.position = "none") +
  scale_y_continuous(position = "right",limits = lims_weight, 
                     breaks = pretty(lims_weight, n = 5)) + labs (y = "")

#Combining everything
boxplot_ant_biometrics_c <- (boxplot_ant_length_c + boxplot_ant_width_c + 
                             boxplot_ant_weight_c) +
  theme( plot.margin      = margin(-10, 20, 5, 20),
         panel.spacing    = unit(1.2, "lines"),
         panel.spacing.y =  unit(0.1, "lines"))      # Abstand zwischen Panels


plot_ant_weight_c <- plot_ant_weight + 
  theme(
  legend.position      = c(0.95, 0.80),           # inwards, top right
  legend.justification = c("right", "center"),
  axis.title.x = element_text(margin = margin (t=15))
)
plot_ant_length_c <- plot_ant_length + theme( 
  axis.title.x = element_text(margin = margin (t=15)))
plot_ant_width_c <- plot_ant_width + theme( 
  axis.title.x = element_text(margin = margin (t=15)))

plot_ant_biometrics_c <- (plot_ant_length_c + plot_ant_width_c +
                            plot_ant_weight_c) +
  theme(
    plot.margin      = margin(5, 20, -5, 20),
    panel.spacing    = unit(0.2, "lines")  )
  

#Histograms
top_row <- (plot_ant_length_c + plot_ant_width_c + plot_ant_weight_c) +
  plot_layout(ncol = 3) &
  theme(plot.margin = margin(t = 5, r = 20, b = 0, l = 20))

# Boxplots 
bottom_row <- (boxplot_ant_length_c + plot_spacer() + boxplot_ant_width_c +
                 plot_spacer() + boxplot_ant_weight_c) +
  plot_layout(ncol = 5, widths = c(1.4, 0.0, 1.4, 0, 1.4)) &
  theme(plot.margin = margin(t = -15, r = 20, b = 5, l = 20))  

All_biomet_plot <- top_row / bottom_row
jpeg(file = "All Biometrics Ants_final.jpg",
     width = 60, height = 32, units = "cm", res = 300)
All_biomet_plot
dev.off()



#==============================================================================#
##### 1.2 Allometry between Ants ##### 
###### 1.2.1 Ant length vs. ants width #####

                          # A. LRM - WIDTH vs. LENGTH 
#-----------------------------------------------------------------------------#
#base model no main effects
ant_biomet_width1 <- lm(log10(ant_width) ~ log10(ant_length), data = df)

summary(ant_biomet_width1)
confint(ant_biomet_width1)
bptest(ant_biomet_width1) 
coeftest(ant_biomet_width1, vcov = vcovHC(ant_biomet_width1, type = "HC3"))

#model with only main effects
ant_biomet_width2 <- lm(log10(ant_width) ~ log10(ant_length) + ant_species
                       + forest_type , data = df)

summary(ant_biomet_width2)
confint(ant_biomet_width2)
bptest(ant_biomet_width2) 
coeftest(ant_biomet_width2, vcov = vcovHC(ant_biomet_width2, type = "HC3"))


#model with interacting effect specis an main effects 
ant_biomet_width3 <- lm(log10(ant_width) ~ log10(ant_length) * ant_species
                     + forest_type , data = df)

summary(ant_biomet_width3)
confint(ant_biomet_width3)
bptest(ant_biomet_width3) 
coeftest(ant_biomet_width3, vcov = vcovHC(ant_biomet_width3, type = "HC3"))


AIC(ant_biomet_width1, ant_biomet_width2, ant_biomet_width3)
anova(ant_biomet_width1, ant_biomet_width2, ant_biomet_width3)

jpeg(file = "Residualplots Ant length vs. width.jpg",
     width = 20, height = 18, units = "cm", res = 300)
par(mfrow = c(2, 2))
plot(ant_biomet_width3)
dev.off()

                         # B. Plot - WIDTH vs. LENGTH 
#-----------------------------------------------------------------------------#
#Defining theme for LRM Plots once - reuse it in every plot
lrm_legend_theme <- theme(
  legend.position      = "bottom",         
  legend.justification = "center",
  legend.direction     = "horizontal",
  legend.box           = "horizontal",
  legend.background    = element_rect(fill = "white", colour = NA),
  legend.box.margin    = margin(0, 0, 0, 0),
  legend.margin        = margin(-5, 0, 0, 0),
  legend.key.size      = unit(0.8, "lines"),
  legend.key.spacing   = unit(0.5, "lines"),
  legend.key.spacing.x = unit(0.8, "lines"),
  legend.text          = element_text(size = 18, face = "italic",  hjust = 0.5),
  legend.title         = element_blank(),
  legend.text.align    = 0.5,
  plot.title           = element_text(hjust = 1, size = 18)
)

# Combining most styles used
lrm_style <- list(
  theme_clean(base_size = 18),       # overall text size
  clean_panel_theme,                 # axes + panel cleanup
  lrm_legend_theme,                  # legend at bottom, italic species
  theme( plot.margin  = margin(t = 5.5, r = 20, b = 10, l = 5.5),
         axis.title.x = element_text(hjust = 0.55),
         axis.title.y = element_text(vjust = -0.5)),
  scale_colour_manual(values = ant_colour ) 
  )

ant_length_width <- ggplot(df, aes(x = ant_length,y = ant_width,
                                    colour = ant_species)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(aes(colour = NULL),   # no species-specific regression
              method = "lm",
              se = TRUE,
              linewidth = 1,
              colour = "grey9") +
  scale_x_log10() + 
  scale_y_log10() +
  lrm_style +                       #  universal lrm style defined in 1.2.1.
  labs( x     = "Ant body length [mm]",
        y     = "Ant head width [mm]",
        color = "Species" )


jpeg(filename = "LRM_ant_length_width.jpg", width = 20, 
     height = 18, units = "cm", res = 300)
ant_length_width
dev.off()


                            
            ######### C. Linear Mixed Models  - WILL NOT USE THIS#######            
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
#NEEDS MASSIV BUG FIX
#predict seems to not hav ant_weight and it seems not to be in the data frame.
#Maybe this issue comes from the lmm model and the previous changes in Intercepter 
#and predicter bc I wanted to have length on the X-achse 
#-----------------------------------------------------------------------------#

newdat_biom_width<- expand.grid(
  ant_length  = seq(min(df$ant_length, na.rm = TRUE),
                    max(df$ant_length, na.rm = TRUE),
                    length.out = 200),
  ant_species = levels(df$ant_species),
  forest_type = "primary",   # reference 
  raid_ID     = NA           # population-level prediction
)

# Dummy-Response-colum added
#enables reading of log10(ant_width) but is not used
#newdat_biom_width$ant_width <- 1

# Vorhersage auf log10-Skala (Modellskala)
newdat_biom_width$log_pred_width <- predict( ant_biomet_width_lmm,
  newdata = newdat_biom_width,
  re.form = NA               # Random Effects = 0 (Population Level)
)


#Für ggplot: x und y vorbereiten
newdat_biom_width$log10_length <- log10(newdat_biom_width$ant_length)
newdat_biom_width$log10_width  <- newdat_biom_width$log_pred_width


#Plot
ggplot(df, aes(x = log10(ant_length), y = log10(ant_width),
               color = ant_species)) +
  
  # Rohdatenpunkte
  geom_point(alpha = 0.3, size = 1.3) +
  
  # LMM-Vorhersagekurven
  geom_line( data = newdat_biom_width,aes( x = log10_length, 
                                           y = log10_width,
              color = ant_species),linewidth = 1.4) +

  scale_colour_manual(values = ant_colour) +
  
  # Theme & Labels
  theme_clean(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 14)
  ) +
  
  labs(
    x = "Ant body length (log10 mm)",
    y = "Ant body width (log10 mm)",
    color = "Species",
    title = "LMM-predicted allometric relationship\n(log10 width ~ log10 length × species)"
  )



###### 1.2.2 Ant length vs. ants weight #####


                         # A. LRM - WEIGHT vs. LENGTH
#-----------------------------------------------------------------------------#
#base model no main effects
ant_biomet_weight1 <- lm(log10(ant_weight*1000) ~ log10(ant_length), data = df)

#model with only main effects
ant_biomet_weight2 <- lm(log10(ant_weight*1000) ~ log10(ant_length) + ant_species
                        + forest_type , data = df)

#model with interacting effect specis an main effects 
ant_biomet_weight3 <- lm(log10(ant_weight*1000) ~ log10(ant_length) * ant_species
                        + forest_type , data = df)

summary(ant_biomet_weight3)
confint(ant_biomet_weight3)
AIC(ant_biomet_weight1, ant_biomet_weight2, ant_biomet_weight3)
anova(ant_biomet_weight1, ant_biomet_weight2, ant_biomet_weight3)


                         # B. Plot  WEIGHT vs. LENGTH
#-----------------------------------------------------------------------------#
ant_length_weight <- ggplot(df, aes(x = ant_length,y = ant_weight * 1000,
                                      colour = ant_species)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(aes(colour = NULL),   # no species-specific regression
    method = "lm",
    se = TRUE,
    linewidth = 1,
    colour = "grey9",
    fill   = "grey66",
    alpha  = 0.9  ) +
  scale_x_log10() + 
  scale_y_log10() +
  lrm_style +                       #  universal lrm style defined in 1.2.1.
  labs( x     = "Ant body length [mm]",
        y     = "Ant weight [mg]",
        color = "Species" )


jpeg(filename = "LRM_ant_length_weight.jpg", width = 20, 
     height = 18, units = "cm", res = 300)
ant_length_weight
dev.off()



###### 1.2.3 Ant weight vs. ants width #####

                        # A. LRM - WEIGHT vs. WIDTH
#-----------------------------------------------------------------------------#
#base model no main effects
ant_biomet_ww1 <- lm(log10(ant_weight*1000) ~ log10(ant_width), data = df)

#model with only main effects
ant_biomet_ww2 <- lm(log10(ant_weight*1000) ~ log10(ant_width) + ant_species
                         + forest_type , data = df)

#model with interacting effect specis an main effects 
ant_biomet_ww3 <- lm(log10(ant_weight*1000) ~ log10(ant_width) * ant_species
                         + forest_type , data = df)

summary(ant_biomet_ww3)
confint(ant_biomet_ww3)
AIC(ant_biomet_ww1, ant_biomet_ww2, ant_biomet_ww3)
anova(ant_biomet_ww1, ant_biomet_ww2, ant_biomet_ww3)


                       # B. Plot  WEIGHT vs. WIDTH
#-----------------------------------------------------------------------------#

ant_width_weight <- ggplot(df, aes(x = ant_width,y = ant_weight * 1000,
                                    colour = ant_species)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(aes(colour = NULL),   # no species-specific regression
              method = "lm",
              se = TRUE,
              linewidth = 1,
              colour = "grey9") +
  scale_x_log10() + 
  scale_y_log10() +
  lrm_style +                       #  universal lrm style defined in 1.2.1.
  labs( x     = "Ant head width [mm]",
        y     = "Ant weight [mg]",
        color = "Species" )


jpeg(filename = "LRM_ant_weight_width.jpg", width = 20, 
     height = 18, units = "cm", res = 300)
ant_width_weight
dev.off()


###### 1.2.4 Combining all LRM Plots #####

ant_length_width_c <- ant_length_width + theme(legend.position = "none")
ant_width_weight_c <- ant_width_weight + theme(legend.position = "none")

lrm_ants <- (ant_length_width_c + ant_length_weight + ant_width_weight_c) +
  boxp_panel

jpeg(filename = "LRM_Panel.jpg", width = 60, 
     height = 18, units = "cm", res = 300)
lrm_ants
dev.off()

#Plot used in Master thesis
Bimetric_ant_resutls_ma <- All_biomet_plot/
  ((ant_length_width_c | plot_spacer()) + plot_layout(widths = c(1, 2)) +
     theme(plot.margin = margin(t = 50))) 

jpeg(file = "Biometrics Ants_MA.jpg",
     width = 60, height = 60, units = "cm", res = 300)
Bimetric_ant_resutls_ma
dev.off()



#######==============================.===================================#######

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

#carrier factor
df_dis <- df_dis %>%
  mutate(carrier_type = if_else(ant_number > 1, "multiple", "single"))


# Make sure its a factor
df_dis$dismembered <- factor(df_dis$dismembered, levels = c("no", "yes"))
df_dis$ant_species <- factor(df_dis$ant_species)
df_dis$forest_type <- factor(df_dis$forest_type)

# Add important factors 
df_dis["prey_shape"] <- df_dis$prey_length / df_dis$prey_width
df_dis["prey_area"] <- df_dis$prey_length * df_dis$prey_width


#=============================================================================#
##### 2.1 Size differences in Prey Items ##### 
#=============================================================================#
           ######### A. Welch t-test ########
#Welch t-test - against FOREST TYP
t1 <- t.test(log10(prey_length) ~ forest_type, data = df_dis)
t2 <- t.test(log10(prey_width) ~ forest_type, data = df_dis)
t3 <- t.test(log10(prey_weight*1000) ~ forest_type, data = df_dis)
t4 <- t.test(log10(prey_area) ~ forest_type, data = df_dis)
t5 <- t.test(log10(prey_shape) ~ forest_type, data = df_dis)

#Welch t-test - against ANT SPECIES
t6  <- t.test(log10(prey_length) ~ ant_species, data = df_dis)
t7  <- t.test(log10(prey_width) ~ ant_species, data = df_dis)
t8  <- t.test(log10(prey_weight*1000) ~ ant_species, data = df_dis)
t9  <- t.test(log10(prey_area) ~ ant_species, data = df_dis)
t10 <- t.test(log10(prey_shape) ~ ant_species, data = df_dis)

#Extracting Results
get_t_results <- function(test, metric, groupvar, data) {
  
  # Get actual factor levels from df
  groups <- levels(as.factor(data[[groupvar]]))
  
  # p-value and stars
  p <- test$p.value
  sig <- ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*",
                              ifelse(p < 0.1, ".", ""))))
  
  data.frame(
    Metric             = metric,
    Grouping           = groupvar,              # use actual df column name
    `Group 1 name`     = groups[1],
    `Group 2 name`     = groups[2],
    `Mean (Group 1)`   = round(test$estimate[1], 3),
    `Mean (Group 2)`   = round(test$estimate[2], 3),
    t                  = round(test$statistic, 3),
    df                 = round(test$parameter, 1),
    p                  = signif(p, 3),
    Significance       = sig,
    Scale              = "log10-transformed"
  )
}

#Build one consitent table
results <- dplyr::bind_rows(
  get_t_results(t1,  "prey length [mm]", "forest_type", df_dis),
  get_t_results(t2,  "prey width [mm]",  "forest_type", df_dis),
  get_t_results(t3,  "prey weight [mg]", "forest_type", df_dis),
  get_t_results(t4,  "prey area [mm²]",  "forest_type", df_dis),
  get_t_results(t5,  "prey shape [-]",   "forest_type", df_dis),
  
  get_t_results(t6,  "prey length [mm]", "ant_species", df_dis),
  get_t_results(t7,  "prey width [mm]",  "ant_species", df_dis),
  get_t_results(t8,  "prey weight [mg]", "ant_species", df_dis),
  get_t_results(t9,  "prey area [mm²]",  "ant_species", df_dis),
  get_t_results(t10, "prey shape [-]",   "ant_species", df_dis)
)
write.csv(results, "prey_ttests_results.csv", row.names = FALSE)


#Checking data with histogarms and Boxplots
par(mfrow = c(2, 3))
par(mar = c(4, 4, 2, 1))
hist(df_dis$prey_length, main = "Prey length", col = "grey",
     xlab = "length [mm]", ylim = c(0, 1000))
hist(df_dis$prey_width, main = "Prey width [mm]",col = "grey",
     xlab = "width [mm]", ylim = c(0, 1000))
hist((1000*df_dis$prey_weight), main = "Prey weight", col = "grey",
     xlab = "Weight [mg]", ylim = c(0, 2500))
hist(df_dis$prey_area, main = expression("Prey area ["*mm^2*+"]"),col = "grey",
     xlab = "width [mm]", ylim = c(0, 1000))
hist(df_dis$prey_shape, main = "Prey shape",col = "grey",
     xlab = "Shape Index (length/width)", ylim = c(0, 3000))


           ######### B. Varianz test ########
# Levene-Test and KS Brown–Forsythe for Varianz tests

#Defining Vektor for Variables
vars_p <- c("prey_length", "prey_width", "prey_weight",
                          "prey_area", "prey_shape")

#Setting up funktion
get_var_tests <- function(var, group, data) {
  # y-Vektor: log10-transformiert, Sonderfall prey_weight * 1000
  if (var == "prey_weight") { y <- log10(data[[var]] * 1000)} 
    else { y <- log10(data[[var]]) }
  
  g <- data[[group]]
  
  # Levene (mean-abs-deviation)
  L  <- leveneTest(y, g)                     
  # Brown–Forsythe (median-abs-deviation)
  BF <- leveneTest(y, g, center = median)    
  
  out <- data.frame(
    metric    = var,
    grouping  = group,
    levene_F  = L$`F value`[1],
    levene_p  = L$`Pr(>F)`[1],
    bf_F      = BF$`F value`[1],
    bf_p      = BF$`Pr(>F)`[1]
  )
  return(out)
}
  
forest_var <- do.call(rbind, lapply(vars_p, get_var_tests, group = "forest_type",
                                      data = df_dis))
species_var <- do.call(rbind, lapply(vars_p, get_var_tests, group = "ant_species",
                                      data = df_dis))
all_var <- bind_rows(forest_var, species_var)
#Fuck Excel...
all_var <- all_var %>% mutate(across(where(is.numeric), ~ round(., 6)))
write.table(all_var, file = "variance_tests_prey_biometrics.csv",
            row.names = FALSE,)



           ######### C. Boxplot Style and Theme #######
# numeric positioning of forest (1 = primary, 2 = secondary)
# with ant specivic scaling 
offset <- 0.22   # settin distance of boxes within forest 

#adding numeric positioning 
df_dis <- df_dis |>mutate(
    forest_center = if_else(forest_type == "primary", 1, 2),      
    species_num = as.numeric(ant_species),          
    x_pos       = forest_center + if_else(species_num == 1,
                                       -offset, offset))

#Defining Style
boxp_style_prey <- list(
  geom_boxplot( aes( x = x_pos, group = interaction(forest_type, ant_species)),
    width          = 0.3,      
    size           = 0.6,      
    colour         = "grey11", # lines box
    alpha          = 0.55,     
    outliers       = FALSE,    
    outlier.shape  = 21,
    outlier.colour = "grey11",
    outlier.fill   = NA,
    outlier.size   = 1,
    outlier.stroke = 0.3,
    outlier.alpha  = 0.3,
    staplewidth    = 0.5,
    show.legend = FALSE,
    position       = position_dodge2(
      width    =  0.65 ,         # distance species boxes within forest typ
      preserve = "single")),
  dorylus_colour,
  theme_clean(base_size = 18) %+replace% theme(
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background   = element_blank(),
    panel.border       = element_blank(),
    plot.title         = element_text(hjust = 1)
  ))

# Pointcloud with Outliners
boxp_distribution_prey <- list(
    geom_sina(
      aes(x = x_pos,             # uses defined position in regard of forest
          group = interaction(forest_type, ant_species)),
      size     = 0.8,             
      alpha    = 0.25,             
      color    = "grey75",         
      method   = "density",
      maxwidth = 0.40              # for point cloud 
    ))

boxp_legend_theme_prey <- list (theme(
    legend.position      = c(0.98, 0.98),           # inwards, top right
    legend.justification = c("right", "top"),       # anker point
    legend.background    = element_rect(fill = "white", colour = NA),
    legend.box.margin    = margin(0, 0, 0, 0),
    legend.margin        = margin(0, 0, 0, 0),
    legend.key.size      = unit(0.7, "lines"),      # compact legend
    legend.key.spacing   = unit(0.5, "lines"),
    legend.key.spacing.x = unit(1.2, "lines"),
    legend.text          = element_text(size = 18, face = "italic"),
    legend.title         =  NULL, 
    legend.text.align  = 0.5,
    legend.title.align = 0.5), 
  guides (colour = "none",        # fill legend only
    fill = guide_legend(
    override.aes = list(
      shape    = 22,
      alpha    = 0.55,
      size     = 5,    
      linetype = 0  )   # no border 
    )))

           ######### D. Generating Boxplots #######

# PREY LENGTH
boxplot_prey_length <- ggplot(df_dis, aes(y = prey_length, fill = ant_species,
                                colour = ant_species)) +
  boxp_distribution_prey +
  boxp_style_prey +
  scale_x_continuous(breaks = c(1, 2),                       # Forest-center
                     labels = c("primary forest", "secondary forest"),
                     name   = NULL) + labs(y = "Prey length [mm]") +
  labs(fill = NULL, colour = NULL) +
  coord_cartesian(ylim = c(0, 20)) +
  boxp_legend_theme_prey            
boxplot_prey_length


# PREY WIDTH
boxplot_prey_width <- ggplot(df_dis, aes(y = prey_width, fill = ant_species,
                                         colour = ant_species)) +
  boxp_distribution_prey +
  boxp_style_prey +
  scale_x_continuous(breaks = c(1, 2),
                     labels = c("primary forest", "secondary forest"),
                     name   = NULL) +
  labs(y = "Prey width [mm]") +
  labs(fill = NULL, colour = NULL) +
  boxp_legend_theme_prey
boxplot_prey_width


# PREY WEIGHT
boxplot_prey_weight <- ggplot(df_dis, aes(y = (prey_weight*1000), fill = ant_species,
                                          colour = ant_species)) +
  boxp_distribution_prey +
  boxp_style_prey +
  scale_x_continuous(breaks = c(1, 2),
                     labels = c("primary forest", "secondary forest"),
                     name   = NULL) +
  labs(y = "Prey weight [mg]") +
  labs(fill = NULL, colour = NULL) +
  coord_cartesian(ylim = c(0, 80)) +
  boxp_legend_theme_prey
boxplot_prey_weight


# PREY AREA
boxplot_prey_area <- ggplot(df_dis, aes(y = prey_area, fill = ant_species,
                                        colour = ant_species)) +
  boxp_distribution_prey +
  boxp_style_prey +
  scale_x_continuous(breaks = c(1, 2),
                     labels = c("primary forest", "secondary forest"),
                     name   = NULL) +
  labs(y = "Prey area [mm²]") +
  labs(fill = NULL, colour = NULL) +
  coord_cartesian(ylim = c(0, 80)) +
  boxp_legend_theme_prey
boxplot_prey_area


# PREY SHAPE
boxplot_prey_shape <- ggplot(df_dis, aes(y = prey_shape, fill = ant_species,
                                         colour = ant_species)) +
  boxp_distribution_prey +
  boxp_style_prey +
  scale_x_continuous(breaks = c(1, 2),
                     labels = c("primary forest", "secondary forest"),
                     name   = NULL) +
  labs(y = "Prey shape index") +
  labs(fill = NULL, colour = NULL) +
  coord_cartesian(ylim = c(0, 10)) +
  boxp_legend_theme_prey
boxplot_prey_shape

           ######### E. Combining Plots #######
#Preperation 
boxplot_prey_length_c <- boxplot_prey_length + boxp_panel +
  theme(legend.position = "none") 
boxplot_prey_width_c  <- boxplot_prey_width + boxp_panel +
  theme(legend.position = "none")
boxplot_prey_weight_c <- boxplot_prey_weight + boxp_panel 
boxplot_prey_area_c   <- boxplot_prey_area + boxp_panel +
  theme(legend.position = "none")
boxplot_prey_shape_c  <- boxplot_prey_shape + boxp_panel +
  theme(legend.position = "none")


#backup - this workded
boxp_prey_combined <- 
  (boxplot_prey_length_c + boxplot_prey_width_c + boxplot_prey_weight_c) /
  (boxplot_prey_area_c + boxplot_prey_shape_c)  + 
  boxp_legend_theme_prey

jpeg(file = "Boxplot Prey All.jpg",
     width = 60, height = 36, units = "cm", res = 300)
boxp_prey_combined
dev.off()

boxp_prey_combined2 <- 
  (boxplot_prey_length_c + boxplot_prey_width_c + boxplot_prey_weight_c) +
  boxp_legend_theme_prey

jpeg(file = "Boxplot Prey All_MA.jpg",
     width = 60, height = 20, units = "cm", res = 300)
boxp_prey_combined2
dev.off()

##### 2.2 Allometry in Prey Items ####
###### 2.2.1 Prey length vs. prey weight ######

                    # A. LRM - LENGTH vs. WEIGHT - PREY
#-----------------------------------------------------------------------------#
#base model no main effects
prey_biomet_weight1 <-lm(log10(prey_weight*1000) ~ log10(prey_length), 
                         data = df_dis)

#model with only main effects
prey_biomet_weight2 <- lm(log10(prey_weight*1000) ~ log10(prey_length) + 
                            ant_species + forest_type , data = df_dis)

#model with interacting effect specis an main effects 
prey_biomet_weight3 <- lm(log10(prey_weight*1000) ~ log10(prey_length) * ant_species
                         + forest_type , data = df_dis)

summary(prey_biomet_weight3)
confint(prey_biomet_weight3)
AIC(prey_biomet_weight1, prey_biomet_weight2, prey_biomet_weight3)
anova(prey_biomet_weight1, prey_biomet_weight2, prey_biomet_weight3)


# B. Plot  WEIGHT vs. LENGTH
#-----------------------------------------------------------------------------#
prey_length_weight <- ggplot(df_dis, aes(x = prey_length,y = prey_weight * 1000,
                                    colour = ant_species)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(aes(colour = NULL),   # no species-specific regression
              method = "lm",
              se = TRUE,
              linewidth = 1,
              colour = "grey9",
              fill   = "grey66",
              alpha  = 0.9  ) +
  scale_x_log10() + 
  scale_y_log10() +
  lrm_style +                       #  universal lrm style defined in 1.2.1.
  labs( x     = "Prey length [mm]",
        y     = "Prey weight [mg]",
        color = "Species" )


jpeg(filename = "LRM_Prey_length_weight.jpg", width = 20, 
     height = 18, units = "cm", res = 300)
prey_length_weight
dev.off()





###### 2.2.2 Prey length vs. prey width ######

                    # A. LRM - LENGTH vs. WIDTH - PREY
#-----------------------------------------------------------------------------#
#base model no main effects
prey_biomet_width1 <-lm(log10(prey_width) ~ log10(prey_length), 
                         data = df_dis)

#model with only main effects
prey_biomet_width2 <- lm(log10(prey_width) ~ log10(prey_length) + 
                            ant_species + forest_type , data = df_dis)

#model with interacting effect specis an main effects 
prey_biomet_width3 <- lm(log10(prey_width) ~ log10(prey_length) * ant_species
                         + forest_type , data = df_dis)

summary(prey_biomet_width3)
confint(prey_biomet_width3)
AIC(prey_biomet_width1, prey_biomet_width2, prey_biomet_width3)
anova(prey_biomet_width1, prey_biomet_width2, prey_biomet_width3)


# B. Plot  WEIGHT vs. LENGTH
#-----------------------------------------------------------------------------#
prey_length_width <- ggplot(df_dis, aes(x = prey_length,y = prey_width,
                                    colour = ant_species)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(aes(colour = NULL),   # no species-specific regression
              method = "lm",
              se = TRUE,
              linewidth = 1,
              colour = "grey9",
              fill   = "grey66",
              alpha  = 0.9  ) +
  scale_x_log10() + 
  scale_y_log10() +
  lrm_style +                       #  universal lrm style defined in 1.2.1.
  labs( x     = "Prey Length [mm]",
        y     = "Prey weight [mg]",
        color = "Species" )


jpeg(filename = "LRM_Prey_length_width.jpg", width = 20, 
     height = 18, units = "cm", res = 300)
prey_length_width
dev.off()





















###### 2.2.3 Prey weight vs. prey width ######


                        # A. LRM - WEIGHT vs. WIDTH
#-----------------------------------------------------------------------------#
#base model no main effects
prey_biomet_ww1 <- lm(log10(prey_weight*1000) ~ log10(prey_width), data = df_dis)

#model with only main effects
prey_biomet_ww2 <- lm(log10(prey_weight*1000) ~ log10(prey_width) + ant_species
                     + forest_type , data = df_dis)

#model with interacting effect specis an main effects 
prey_biomet_ww3 <- lm(log10(prey_weight*1000) ~ log10(prey_width) * ant_species
                     + forest_type , data = df_dis)

summary(prey_biomet_ww3)
confint(prey_biomet_ww3)
AIC(prey_biomet_ww1, prey_biomet_ww2, prey_biomet_ww3)
anova(prey_biomet_ww1, prey_biomet_ww2, prey_biomet_ww3)


                     # B. Plot  WEIGHT vs. WIDTH -Prey
#-----------------------------------------------------------------------------#

prey_width_weight <- ggplot(df_dis, aes(x = prey_width,y = prey_weight * 1000,
                                   colour = ant_species)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(aes(colour = NULL),   # no species-specific regression
              method = "lm",
              se = TRUE,
              linewidth = 1,
              colour = "grey9") +
  scale_x_log10() + 
  scale_y_log10() +
  lrm_style +                       #  universal lrm style defined in 1.2.1.
  labs( x     = "Prey width [mm]",
        y     = "Prey weight [mg]",
        color = "Species" )


jpeg(filename = "LRM_Prey_weight_width.jpg", width = 20, 
     height = 18, units = "cm", res = 300)
prey_width_weight
dev.off()

###### 2.2.4 Prey weight vs. prey area ######


                      # A. LRM - WEIGHT vs. AREA
#-----------------------------------------------------------------------------#
#base model no main effects
prey_biomet_dens1 <- lm(log10(prey_weight*1000) ~ log10(prey_area), data = df_dis)

#model with only main effects
prey_biomet_dens2 <- lm(log10(prey_weight*1000) ~ log10(prey_area) + ant_species
                      + forest_type , data = df_dis)

#model with interacting effect specis an main effects 
prey_biomet_dens3 <- lm(log10(prey_weight*1000) ~ log10(prey_area) * ant_species
                      + forest_type , data = df_dis)

summary(prey_biomet_dens3)
confint(prey_biomet_dens3)
AIC(prey_biomet_dens1, prey_biomet_dens2, prey_biomet_dens3)
anova(prey_biomet_dens1, prey_biomet_dens2, prey_biomet_dens3)


# B. Plot  WEIGHT vs. WIDTH -Prey
#-----------------------------------------------------------------------------#

prey_density <- ggplot(df_dis, aes(x = prey_area,y = prey_weight * 1000,
                                    colour = ant_species)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(aes(colour = NULL),   # no species-specific regression
              method = "lm",
              se = TRUE,
              linewidth = 1,
              colour = "grey9") +
  scale_x_log10() + 
  scale_y_log10() +
  lrm_style +                       #  universal lrm style defined in 1.2.1.
  labs( x     = "Prey area [mm²]",
        y     = "Prey weight [mg]",
        color = "Species" )


jpeg(filename = "LRM_Prey_weight_area.jpg", width = 20, 
     height = 18, units = "cm", res = 300)
prey_density
dev.off()

#==================================.===========================================

                   #### 3 Ant-Prey-Relationship  #### 

#=============================================================================#

###### 3.1 Ant weight vs prey weight (single workers) ######
#=============================================================================#
####### Data preparation #####
#setting species and forest as factor (aka categorical variable)
df_single$ant_species <- factor(df_single$ant_species)
df_single$forest_type <- factor(df_single$forest_type)

#In prey_shape and ant_size are INF (when not log transformed data == 1.00)
#We assume the index is just 0 here to get rid of errors 
df_single <- df_single %>%
  mutate(across(c(ant_size, prey_shape),
                ~ ifelse(is.infinite(.), 0, .))) #if Inf/-Inf write 0 

####### 3.1.1 LRM comparison  #####
#base model without explaining variables 
#weight relations for single carriers Nr.1
weight_sc1 <- lm(log10(prey_weight*1000) ~ log10(ant_weight*1000),
                 data = df_single) 

summary(weight_sc1 )
confint(weight_sc1 ) # 95%-Convidence intervall
bptest(weight_sc1) 
coeftest(weight_sc1, vcov = vcovHC(weight_sc1, type = "HC3"))


#model with only main effects
weight_sc2 <- lm(log10(prey_weight*1000) ~ log10(ant_weight*1000) + 
                   ant_species + forest_type +  prey_shape, data = df_single)

summary(weight_sc2)
confint(weight_sc2)

#AIC(weight_sc1 , weight_sc2)  #compares models - Akaike Information Criterion
#plot(weight_sc2) #standard plots to check homoskedasticity, linearity, leverage

#Breusch-Pagan-test on homoskedasticity (is Varianz in Residuals constant)
bptest(weight_sc2) 
coeftest(weight_sc2, vcov = vcovHC(weight_sc2, type = "HC3"))


#model with interacting effect specis an main effects 
#test this to see if slope is different when species are considered
weight_sc3 <- lm(log10(prey_weight*1000) ~ log10(ant_weight*1000) * ant_species+ 
                   forest_type +  prey_shape, data = df_single)

summary(weight_sc3)
confint(weight_sc3)
bptest(weight_sc3) 
coeftest(weight_sc3, vcov = vcovHC(weight_sc3, type = "HC3"))

AIC(weight_sc1, weight_sc2, weight_sc3)
anova(weight_sc1, weight_sc2, weight_sc3)



####### 3.1.2 LMM calculation #####
#-----------------------------------------------------------------------------#
#Creating a LMM based on LRM with only main/fixed effects 
#This includes possible random variation due to the raids (adds random effect)
#We do this based on the rusltus of comparison between the models above


#The random intercept for raid identity explained 13% of the variance (ICC=0.13)

weight_sc_LMM <- lmer(log10(prey_weight*1000) ~ log10(ant_weight*1000) + 
                        ant_species + forest_type + prey_shape + (1 | raid_ID),
                      data = df_single)

summary(weight_sc_LMM)
AIC(weight_sc1, weight_sc2, weight_sc3, weight_sc_LMM)

#checking variation within Raids
VarCorr(weight_sc_LMM)
icc(weight_sc_LMM)

                              #  DHARMa 
#-----------------------------------------------------------------------------#

sim_weight_sc_LMM <- simulateResiduals(weight_sc_LMM)  #Simulating Residuals

# plot(sim_weight_sc_LMM ) #plot diagnostic

# Testing Distribution & Dispersion
testDispersion(sim_weight_sc_LMM )
testUniformity(sim_weight_sc_LMM )
testResiduals(sim_weight_sc_LMM )  



####### 3.1.3 LMM Plot #####
#-----------------------------------------------------------------------------#
# Sequencing Ant weight to have fixed points for the predict grid
# This is done in the original scale and later transformed from g to mg
ant_seq_weight_sc_LMM <- seq(
  from = min((df_single$ant_weight), na.rm = TRUE),
  to   = max((df_single$ant_weight), na.rm = TRUE),
  length.out = 100
)

#Prediction-Grid (for both spices, primary forest, mean prey_shape)
#To calculate predicts and generate one line all effects net to be standardized
newdat_weight_sc_LM <- expand.grid(
  ant_weight   = ant_seq_weight_sc_LMM,
  ant_species  = levels(df_single$ant_species),
  forest_type  = "primary",
  prey_shape   = mean(df_single$prey_shape, na.rm = TRUE)
)

#LMM-Predicts (only Fixed Effects)
newdat_weight_sc_LM$pred_log_weight_sc_LMM <- predict(
  weight_sc_LMM,
  newdata = newdat_weight_sc_LM,
  re.form = NA   # no random effects in line
)

# Gemeinsames Legendendesign für LRM + LMM Plots
lmm_legend_theme <- list(dorylus_colour,
  theme(
    legend.position      = c(0.98, 0.07),   # im Panel: rechts, über x-Achse
    legend.justification = c(1, 0),         # unten-rechte Ecke "ankern"
    legend.direction     = "vertical",      # Einträge untereinander
    legend.box           = "vertical",
    legend.background    = element_rect(fill = "white", colour = NA),
    legend.key           = element_rect(fill = NA, colour = NA),
    legend.key.width     = unit(0.4, "cm"),
    legend.key.height    = unit(0.15, "cm"),
    legend.key.size      = unit(0.25, "lines"),
    legend.text          = element_text(size = 14, face = "italic", hjust = 0),
    legend.title         = element_blank(),
    legend.text.align    = 0,
    legend.box.margin    = margin(0, 0, 0, 0),
    legend.margin        = margin(0, 0, 0, 0)),
  guides(
    colour = guide_legend(
      override.aes = list(
        shape    = 21,     # gefülltes Quadrat
        size     = 2,    # Icon-Größe
        linetype = 0,       # keine Linie, nur Symbol
        alpha  = 0.7))))


plot_weight_sc_LMM <- ggplot(df_single, aes(x = ant_weight * 1000,
      y = prey_weight * 1000, colour = ant_species, fill = ant_species)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_line( data = newdat_weight_sc_LMM,
    aes( x = ant_weight * 1000,
         y = 10^pred_log_weight_sc_LMM,   # log10(mg) -> mg
      colour = ant_species), linewidth = 1) +
  
  scale_x_log10(name = "Ant weight [mg]") +
  scale_y_log10(name = "Prey weight [mg]") +
  labs(colour = "Dorylus",fill   = "Dorylus") +
  lrm_style +
  lmm_legend_theme

plot_weight_sc_LMM

# file save
jpeg(file = "Ant vs Prey Weight LMM.jpg",
     width = 25, height = 18, units = "cm", res = 300)
plot_weight_sc_LMM
dev.off()


#=============================================================================#


####### 3.1.4 LRM Plot - most basic #####

plot_weight_sc_basic <- ggplot(df_single, aes(x = log10(1000*ant_weight),
                                              y = log10(1000*prey_weight),
                                              color = ant_species)) +
  #sets points  
  geom_point(alpha = 0.7, size = 0.7) +  
  
  #regression line, se sets confidence area 
  geom_smooth (method = "lm", se = TRUE, linewidth = .7)+  
  
  #add horizontel reference line   
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, color = "red") +
  labs(
    x = "Ant weight (log10, mg)",
    y = "Prey weight (log10, mg)",
    color = "Dorylus"
  ) +
  dorylus_colour +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)) 

# file save
jpeg(file = "Ant vs Prey Weight LRM Basic.jpg",
     width = 20, height = 18, units = "cm", res = 300)
plot_weight_sc_basic
dev.off()
















####### 3.1.5 LRM Plot - for best fitting model (weight_sc2) ######
        #------------Lrgend needs fixing in thi----------#
# LRM model 3 prediction
pred_weight_LRM <- predict(weight_sc2, newdata  = newdat_weight_sc_LM,
                   interval = "confidence")

newdat_weight_sc_LM$fit_log_prey <- pred_weight_LRM[, "fit"]
newdat_weight_sc_LM$lwr_log_prey <- pred_weight_LRM[, "lwr"]
newdat_weight_sc_LM$upr_log_prey <- pred_weight_LRM[, "upr"]


#confidanc of regressionlines
rim_LM_weight <-   geom_ribbon(data = newdat_weight_sc_LM,aes(  
  x = ant_weight * 1000,
  ymin = 10^lwr_log_prey,
  ymax = 10^upr_log_prey,
  fill = ant_species ),
  alpha = 0.4,
  colour = NA,
  inherit.aes = FALSE)

plot_weight_sc_LM <- ggplot(
  df_single, aes( x = ant_weight * 1000,   #log transformation later in code
                  y = prey_weight * 1000, colour = ant_species, 
                  fill   = ant_species)) +
  geom_point(alpha = 0.3, size = 1.5) +    #adds data points    
  
  geom_line( data = newdat_weight_sc_LM,          # Regressionline (LRM Model 2)
    aes(x = ant_weight * 1000,y = 10^fit_log_prey # zurücktransformiert in mg
    ), linewidth = 1) +
  
  #rim_LM_weight +

  scale_x_log10(name = "Ant weight [mg]") +
  scale_y_log10(name = "Prey weight [mg]")+
  labs(colour = "Dorylus", fill   = "Dorylus")+
  lrm_style +          # defined in chapter 1
  lmm_legend_theme +# defined under 3.1 LMM
  guides( colour = guide_legend(          # gleiche Icons wie im LMM
      override.aes = list(
        shape    = 21,   # Punkte
        size     = 3.5)))
plot_weight_sc_LM 

# file save
jpeg(file = "Ant vs Prey Weight LRM.jpg",
     width = 20, height = 18, units = "cm", res = 300)
plot_weight_sc_LM
dev.off()

jpeg(file = "Residualplots Ant vs Prey Weight.jpg",
     width = 20, height = 18, units = "cm", res = 300)
par(mfrow = c(2, 2))
plot(weight_sc2)
dev.off()

##### 3.2 Dimensional matching  #### 
#=============================================================================#
###### 3.2.1 Prey area  vs. Ant area ###### 
#=============================================================================#
                ####### A LRM comparison – Area  #####
df_single["prey_area"] <- df_single$prey_length * df_single$prey_width
df_single["ant_size"] <- df_single$ant_length * df_single$ant_width

# base model without explanatory variables
area_sc1 <- lm(log10(prey_area) ~ log10(ant_size),
               data = df_single)

summary(area_sc1)
confint(area_sc1)
bptest(area_sc1)

# model with only main effects
area_sc2 <- lm(log10(prey_area) ~ log10(ant_size) +
                 ant_species + forest_type + prey_shape,
               data = df_single)

summary(area_sc2)
confint(area_sc2)
#AIC(area_sc1, area_sc2)
bptest(area_sc2)


# interaction model (species × ant_size)
area_sc3 <- lm(log10(prey_area) ~ log10(ant_size) * ant_species +
                 forest_type + prey_shape,
               data = df_single)

summary(area_sc3)
confint(area_sc3)
bptest(area_sc3)
AIC(area_sc1, area_sc2, area_sc3)
anova(area_sc1, area_sc2, area_sc3)



#=============================================================================#
                ####### B. LMM calculation – Area #####

area_sc_LMM <- lmer(log10(prey_area) ~ log10(ant_size) +
                      ant_species + forest_type + prey_shape +
                      (1 | raid_ID),
                    data = df_single)

summary(area_sc_LMM)
AIC(area_sc1, area_sc2, area_sc3, area_sc_LMM)

VarCorr(area_sc_LMM)
icc(area_sc_LMM)


# DHARMa diagnostics
sim_area_sc_LMM <- simulateResiduals(area_sc_LMM)
#plot(sim_area_sc_LMM)
testDispersion(sim_area_sc_LMM)
testUniformity(sim_area_sc_LMM)
testResiduals(sim_area_sc_LMM)



#=============================================================================#
                ####### C. LMM Plot – Ant area vs Prey area #####

# sequence for predictions
ant_seq_area_sc_LMM <- seq(
  from = min(df_single$ant_size, na.rm = TRUE),
  to   = max(df_single$ant_size, na.rm = TRUE),
  length.out = 100
)

# prediction grid
newdat_area_sc_LM <- expand.grid(
  ant_size    = ant_seq_area_sc_LMM,
  ant_species = levels(df_single$ant_species),
  forest_type = "primary",
  prey_shape  = mean(df_single$prey_shape, na.rm = TRUE)
)

# fixed-effect predictions
newdat_area_sc_LM$pred_log_area_sc_LMM <- predict(
  area_sc_LMM,
  newdata = newdat_area_sc_LM,
  re.form = NA
)

plot_area_sc_LMM <- ggplot(
  df_single,
  aes(x = ant_size,
      y = prey_area,
      colour = ant_species)
) +
  geom_point(alpha = 0.7, size = 0.7) +
  geom_line(
    data = newdat_area_sc_LM,
    aes(
      x = ant_size,
      y = 10^pred_log_area_sc_LMM,
      colour = ant_species
    ),
    linewidth = 1
  ) +
  scale_x_log10(name = "Ant area [mm²]") +
  scale_y_log10(name = "Prey area [mm²]") +
  labs(colour = "Dorylus", fill = "Dorylus") +
  lrm_style +
  lmm_legend_theme

plot_area_sc_LMM

jpeg(file = "Ant vs Prey Area LMM.jpg",
     width = 25, height = 18, units = "cm", res = 300)
plot_area_sc_LMM
dev.off()



#=============================================================================#
                ####### D. LRM Plot – best fitting model (area_sc3) #####

area_sc1 <- lm(log10(prey_area) ~ log10(ant_size), data = df_single)

plot_area_basic_lm <- ggplot(df_single, aes(log10(ant_size), log10(prey_area))) +
  
  geom_point(alpha = .7, size = .7) +
  geom_smooth (method = "lm", se = FALSE, linewidth = 1, color = ant_species)+
  clean_theme +
  labs( x = "Ant area ",
    y = "Prey area)",
    color = "Dorylus") +
  lrm_legend_theme+
  lrm_style
  
plot_area_basic_lm



pred_area_LRM <- predict(
  area_sc3,
  newdata  = newdat_area_sc_LM,
  interval = "confidence"
)

newdat_area_sc_LM$fit_log_prey_area <- pred_area_LRM[, "fit"]
newdat_area_sc_LM$lwr_log_prey_area <- pred_area_LRM[, "lwr"]
newdat_area_sc_LM$upr_log_prey_area <- pred_area_LRM[, "upr"]

plot_area_sc_LM <- ggplot(
  df_single,
  aes(
    x = ant_size,
    y = prey_area,
    colour = ant_species
  )
) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_line(
    data = newdat_area_sc_LM,
    aes(
      x = ant_size,
      y = 10^fit_log_prey_area
    ),
    linewidth = 1
  ) +
  geom_ribbon(
    data = newdat_area_sc_LM,
    aes(
      x = ant_size,
      ymin = 10^lwr_log_prey_area,
      ymax = 10^upr_log_prey_area,
      fill = ant_species
    ),
    alpha = 0.4,
    colour = NA,
    inherit.aes = FALSE
  ) +
  scale_x_log10(name = "Ant area [mm²]") +
  scale_y_log10(name = "Prey area [mm²]") +
  labs(colour = "Dorylus", fill = "Dorylus") +
  lrm_style +
  lmm_legend_theme

plot_area_sc_LM

jpeg(file = "Ant vs Prey Area LRM.jpg",
     width = 20, height = 18, units = "cm", res = 300)
plot_area_sc_LM
dev.off()




###### 3.2.2 Prey shape vs. Ant area ###### 
                ####### A LRM comparison – Shape  #####

# base model without explanatory variables
shape_sc1 <- lm(log10(prey_area) ~ log10(ant_size),
                data = df_single)

summary(shape_sc1)
confint(shape_sc1)
bptest(shape_sc1)

# model with only main effects
shape_sc2 <- lm(log10(prey_area) ~ log10(ant_size) +
                  ant_species + forest_type + prey_shape,
                data = df_single)

summary(shape_sc2)
confint(shape_sc2)
bptest(shape_sc2)


# interaction model (species × ant_size)
shape_sc3 <- lm(log10(prey_area) ~ log10(ant_size) * ant_species +
                  forest_type + prey_shape,
                data = df_single)

summary(shape_sc3)
confint(shape_sc3)
bptest(shape_sc3)
AIC(shape_sc1, shape_sc2, shape_sc3)
anova(shape_sc1, shape_sc2, shape_sc3)



#=============================================================================#
                ####### B. LMM calculation – Shape #####

shape_sc_LMM <- lmer(log10(prey_area) ~ log10(ant_size) +
                       ant_species + forest_type + prey_shape +
                       (1 | raid_ID),
                     data = df_single)

summary(shape_sc_LMM)
AIC(shape_sc1, shape_sc2, shape_sc3, shape_sc_LMM)

VarCorr(shape_sc_LMM)
icc(shape_sc_LMM)


# DHARMa diagnostics
sim_shape_sc_LMM <- simulateResiduals(shape_sc_LMM)
testDispersion(sim_shape_sc_LMM)
testUniformity(sim_shape_sc_LMM)
testResiduals(sim_shape_sc_LMM)



#=============================================================================#
                ####### C. LMM Plot – Ant shape vs Prey shape #####

# sequence for predictions
ant_seq_shape_sc_LMM <- seq(
  from = min(df_single$ant_size, na.rm = TRUE),
  to   = max(df_single$ant_size, na.rm = TRUE),
  length.out = 100
)

# prediction grid
newdat_shape_sc_LM <- expand.grid(
  ant_size    = ant_seq_shape_sc_LMM,
  ant_species = levels(df_single$ant_species),
  forest_type = "primary",
  prey_shape  = mean(df_single$prey_shape, na.rm = TRUE)
)

# fixed-effect predictions
newdat_shape_sc_LM$pred_log_shape_sc_LMM <- predict(
  shape_sc_LMM,
  newdata = newdat_shape_sc_LM,
  re.form = NA
)

plot_shape_sc_LMM <- ggplot(
  df_single,
  aes(x = ant_size,
      y = prey_area,
      colour = ant_species)
) +
  geom_point(alpha = 0.7, size = 0.7) +
  geom_line(
    data = newdat_shape_sc_LM,
    aes(
      x = ant_size,
      y = 10^pred_log_shape_sc_LMM,
      colour = ant_species
    ),
    linewidth = 1
  ) +
  scale_x_log10(name = "Ant shape [mm²]") +
  scale_y_log10(name = "Prey shape [mm²]") +
  labs(colour = "Dorylus", fill = "Dorylus") +
  lrm_style +
  lmm_legend_theme

plot_shape_sc_LMM

jpeg(file = "Ant vs Prey Shape LMM.jpg",
     width = 25, height = 18, units = "cm", res = 300)
plot_shape_sc_LMM
dev.off()



#=============================================================================#
                ####### D. LRM Plot – best fitting model (shape_sc3) #####

shape_sc1 <- lm(log10(prey_area) ~ log10(ant_size), data = df_single)

plot_shape_basic_lm <- ggplot(df_single, aes(log10(ant_size), log10(prey_area), 
                                             color = ant_species)) +
  
  geom_point(alpha = .7, size = .7) +
  geom_smooth(aes(colour = ant_species), method = "lm", se = FALSE, 
              linewidth = 1) +
  clean_theme +
  labs(x = "Ant shape",
       y = "Prey shape",
       color = "Dorylus") +
  lrm_legend_theme +
  lrm_style

plot_shape_basic_lm

jpeg(file = "Ant vs Prey Shape LRM Basic.jpg",
     width = 20, height = 18, units = "cm", res = 300)
plot_shape_sc_LM
dev.off()





pred_shape_LRM <- predict(
  shape_sc3,
  newdata  = newdat_shape_sc_LM,
  interval = "confidence"
)

newdat_shape_sc_LM$fit_log_prey_shape <- pred_shape_LRM[, "fit"]
newdat_shape_sc_LM$lwr_log_prey_shape <- pred_shape_LRM[, "lwr"]
newdat_shape_sc_LM$upr_log_prey_shape <- pred_shape_LRM[, "upr"]

plot_shape_sc_LM <- ggplot(
  df_single,
  aes(
    x = ant_size,
    y = prey_area,
    colour = ant_species
  )
) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_line(
    data = newdat_shape_sc_LM,
    aes(
      x = ant_size,
      y = 10^fit_log_prey_shape
    ),
    linewidth = 1
  ) +
  geom_ribbon(
    data = newdat_shape_sc_LM,
    aes(
      x = ant_size,
      ymin = 10^lwr_log_prey_shape,
      ymax = 10^upr_log_prey_shape,
      fill = ant_species
    ),
    alpha = 0.4,
    colour = NA,
    inherit.aes = FALSE
  ) +
  scale_x_log10(name = "Ant shape [mm²]") +
  scale_y_log10(name = "Prey shape [mm²]") +
  labs(colour = "Dorylus", fill = "Dorylus") +
  lrm_style +
  lmm_legend_theme

plot_shape_sc_LM

jpeg(file = "Ant vs Prey Shape LRM.jpg",
     width = 20, height = 18, units = "cm", res = 300)
plot_shape_sc_LM
dev.off()












#=============================================================================#
##### 3.3 Loading in Single Workers  #### 
#=============================================================================#
###### 3.3.1 Relative load vs. Ant weight (single workers) ######

                ####### A. Linear Regression Models ####
#-----------------------------------------------------------------------------#
#base model no main effects
rel_load1w <- lm(log10(relativ_loading) ~ log10(1000*ant_weight),
                 data = df_single) 

summary(rel_load1w)
confint(rel_load1w) # 95%-Convidence intervall
#plot(rel_load1w) 
bptest(rel_load1w)
coeftest(rel_load1w, vcov = vcovHC(rel_load1w, type = "HC3"))

#model with only main effects
rel_load2w <- lm(log10(relativ_loading) ~ log10(1000*ant_weight) + ant_species+
                   forest_type +  prey_shape, data = df_single)

summary(rel_load2w)
confint(rel_load2w)
#plot(rel_load2w)
bptest(rel_load2w)
coeftest(rel_load2w, vcov = vcovHC(rel_load2w, type = "HC3"))

#model with interacting effect specis an main effects 
rel_load3w <- lm(log10(relativ_loading) ~ log10(1000*ant_weight) * ant_species+ 
                   forest_type + prey_shape, data = df_single)


summary(rel_load3w)
confint(rel_load3w)
bptest(rel_load3w)
coeftest(rel_load3w, vcov = vcovHC(rel_load3w, type = "HC3"))
AIC(rel_load1w, rel_load2w, rel_load3w)
anova(rel_load1w, rel_load2w, rel_load3w)
#plot(rel_load3w)

                ####### B. Linear Mixed Models(not used-probab. need fix)#####
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



                ####### C. Plot for LMM ####
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
jpeg(file = "Relative Load vs. Ant Weight LMM.jpg",
     width = 25, height = 18, units = "cm", res = 300)
plot_rel_weight
dev.off()

                ####### D. Plot for LRM`s ####
#-----------------------------------------------------------------------------#
rel_load1 <- lm(log10(relativ_loading) ~ log10(1000*ant_weight), 
                data = df_single)

# --- Basisplot ---
plot_rel_load_basic_lm <- ggplot(df_single,
  aes(log10(1000*ant_weight), log10(relativ_loading), colour = ant_species)) +
  geom_point(alpha = .7, size = .7) +
  geom_smooth(aes(colour = ant_species), method = "lm", se = FALSE,
              linewidth = 1) +
  clean_theme +
  labs(x = "Ant weight",
       y = "Relative loading",
       colour = "Dorylus") +
  lrm_legend_theme +
  lrm_style

plot_rel_load_basic_lm

jpeg(file = "Relative Load vs. Ant Weight LRM Basic.jpg",
     width = 20, height = 18, units = "cm", res = 300)
plot_rel_load_basic_lm
dev.off()


#  Vorhersagen (for model 2 - fittest ) 
pred_rel_load_LRM <- predict(
  rel_load2w,
  newdata = newdat_rel_loadw,
  interval = "confidence"
)

newdat_rel_loadw$fit_log_rel_load <- pred_rel_load_LRM[, "fit"]
newdat_rel_loadw$lwr_log_rel_load <- pred_rel_load_LRM[, "lwr"]
newdat_rel_loadw$upr_log_rel_load <- pred_rel_load_LRM[, "upr"]


# --- Plot mit Regressionslinien + Bändern ---
plot_rel_load_weightLM <- ggplot( df_single, aes( x = ant_weight,
                                                  y = relativ_loading,
                                                  colour = ant_species)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_line(data = newdat_rel_loadw, aes(x = ant_weight,
                                         y = 10^fit_log_rel_load,
                                         colour = ant_species),
            linewidth = 1,
            inherit.aes = FALSE) +
  geom_ribbon(data = newdat_rel_loadw, aes(
      x = ant_weight,
      ymin = 10^lwr_log_rel_load,
      ymax = 10^upr_log_rel_load,
      fill = ant_species),
    alpha = 0.4,
    colour = NA,
    inherit.aes = FALSE) +
  scale_x_log10(name = "Ant weight") +
  scale_y_log10(name = "Relative loading") +
  labs(colour = "Dorylus", fill = "Dorylus") +
  lrm_style +
  lmm_legend_theme

plot_rel_load_weightLM

jpeg(file = "Relative Load vs. Ant Weight LRM.jpg",
     width = 20, height = 18, units = "cm", res = 300)
plot_rel_load_weightLM
dev.off()

jpeg(file = "Residualplots Relartiv load in single workers.jpg",
     width = 20, height = 18, units = "cm", res = 300)
par(mfrow = c(2, 2))
plot(rel_load2w)
dev.off()


#_____________________________________________________________________________#

###### 3.3.2 Relative load vs. Ant size (single workers) ######
# Relative load is defined by ant_load / ant-weight

                ####### A. Linear Regression Models ####
#-----------------------------------------------------------------------------#
# base model, no main effects
rel_load1s <- lm(log10(relativ_loading) ~ log10(ant_size), data = df_single) 

summary(rel_load1s)
confint(rel_load1s)
bptest(rel_load1s)

# model with only main effects
rel_load2s <- lm(log10(relativ_loading) ~ log10(ant_size) + ant_species +
                   forest_type + prey_shape, data = df_single)

summary(rel_load2s)
confint(rel_load2s)
bptest(rel_load2s)

# model with interaction species × ant_size
rel_load3s <- lm(log10(relativ_loading) ~ log10(ant_size) * ant_species + 
                   forest_type + prey_shape, data = df_single)

summary(rel_load3s)
confint(rel_load3s)
bptest(rel_load3s)
AIC(rel_load1s, rel_load2s, rel_load3s)
anova(rel_load1s, rel_load2s, rel_load3s)


                ####### B. Linear Mixed Model ####
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


                            # DHARMa 
#-----------------------------------------------------------------------------#

sim_rel_load_LMMs <- simulateResiduals(rel_load_LMMs)

# plot(sim_rel_load_LMMs)

testDispersion(sim_rel_load_LMMs)
testUniformity(sim_rel_load_LMMs)
testResiduals(sim_rel_load_LMMs)


                ####### C. Plot for LMM #####
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



                ####### D. Plot for LRM`s ####
#-----------------------------------------------------------------------------#
rel_load1s <- lm(log10(relativ_loading) ~ log10(ant_size), data = df_single)

# --- Basisplot ---
plot_rel_load_basic_lm <- ggplot(
  df_single,
  aes(log10(ant_size), log10(relativ_loading), colour = ant_species)
) +
  geom_point(alpha = .7, size = .7) +
  geom_smooth(aes(colour = ant_species), method = "lm", se = FALSE,
              linewidth = 1) +
  clean_theme +
  labs(
    x = "Ant size",
    y = "Relative loading",
    colour = "Dorylus"
  ) +
  lrm_legend_theme +
  lrm_style

plot_rel_load_basic_lm

jpeg(file = "Relative Load vs. Ant Size LRM Basic.jpg",
     width = 20, height = 18, units = "cm", res = 300)
plot_rel_load_basic_lm
dev.off()


#  Vorhersagen (for model 2 - fittest) 
pred_rel_load_LRM <- predict(
  rel_load2s,
  newdata = newdat_rel_loads,
  interval = "confidence"
)

newdat_rel_loads$fit_log_rel_load <- pred_rel_load_LRM[, "fit"]
newdat_rel_loads$lwr_log_rel_load <- pred_rel_load_LRM[, "lwr"]
newdat_rel_loads$upr_log_rel_load <- pred_rel_load_LRM[, "upr"]


# --- Plot mit Regressionslinien + Bändern ---
plot_rel_load_LM <- ggplot(
  df_single,
  aes(
    x = ant_size,
    y = relativ_loading,
    colour = ant_species
  )
) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_line(
    data = newdat_rel_loads,
    aes(
      x = ant_size,
      y = 10^fit_log_rel_load,
      colour = ant_species
    ),
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  geom_ribbon(
    data = newdat_rel_loads,
    aes(
      x = ant_size,
      ymin = 10^lwr_log_rel_load,
      ymax = 10^upr_log_rel_load,
      fill = ant_species
    ),
    alpha = 0.4,
    colour = NA,
    inherit.aes = FALSE
  ) +
  scale_x_log10(name = "Ant size") +
  scale_y_log10(name = "Relative loading") +
  labs(colour = "Dorylus", fill = "Dorylus") +
  lrm_style +
  lmm_legend_theme

plot_rel_load_LM

jpeg(file = "Relative Load vs. Ant Size LRM.jpg",
     width = 20, height = 18, units = "cm", res = 300)
plot_rel_load_LM
dev.off()

#_____________________________________________________________________________#


###### 3.3.3 Relative load vs. Ant length (single workers) ######
####### A. Linear Regression Models ####
#-----------------------------------------------------------------------------#
# base model – no additional main effects
rel_load1_l <- lm(log10(relativ_loading) ~ log10(ant_length),
                  data = df_single)

summary(rel_load1_l)
confint(rel_load1_l)           # 95% confidence intervals
# plot(rel_load1_l)
bptest(rel_load1_l)
coeftest(rel_load1_l, vcov = vcovHC(rel_load1_l, type = "HC3"))

# model with main effects
rel_load2_l <- lm(log10(relativ_loading) ~ log10(ant_length) +
                    ant_species + forest_type + prey_shape,
                  data = df_single)

summary(rel_load2_l)
confint(rel_load2_l)
# plot(rel_load2_l)
bptest(rel_load2_l)
coeftest(rel_load2_l, vcov = vcovHC(rel_load2_l, type = "HC3"))

# selected model: interaction between ant length and species + main effects 
rel_load3_l <- lm(log10(relativ_loading) ~ log10(ant_length) * ant_species +
                    forest_type + prey_shape,
                  data = df_single)

summary(rel_load3_l)
confint(rel_load3_l)
bptest(rel_load3_l)
coeftest(rel_load3_l, vcov = vcovHC(rel_load3_l, type = "HC3"))

# model comparison (basis für Modellwahl: rel_load3_l)
AIC(rel_load1_l, rel_load2_l, rel_load3_l)
anova(rel_load1_l, rel_load2_l, rel_load3_l)
# plot(rel_load3_l)



####### B. Plots for LRM`s (selected model: rel_load3_l) ####
#-----------------------------------------------------------------------------#
# Basic plot 
plot_rel_load_basic_lm <- ggplot(df_single,aes(x = log10(ant_length),
                                               y = log10(relativ_loading),
                                               colour = ant_species)) +
  geom_point(alpha = 0.7, size = 0.7) +
  geom_smooth(aes(colour = ant_species),
              method = "lm", se = FALSE, linewidth = 1) +
  clean_theme +
  labs(x = "Ant length (log10 mm)",
       y = "Relative loading (log10 prey mass / ant mass)",
       colour = "Dorylus") +
  lrm_legend_theme +
  lrm_style

plot_rel_load_basic_lm

jpeg(file = "Relative Load vs. Ant Length LRM Basic.jpg",
     width = 20, height = 18, units = "cm", res = 300)
plot_rel_load_basic_lm
dev.off()


# Prediction for chosen model rel_load3_l
# Sequenz for Ant length 
ant_seq_rel_load_l <- seq(
  from = min(df_single$ant_length, na.rm = TRUE),
  to   = max(df_single$ant_length, na.rm = TRUE),
  length.out = 100)

# Prediction-Grid (for both species, primary forest, and mean prey_shape)
newdat_rel_load_l <- expand.grid(
  ant_length  = ant_seq_rel_load_l,
  ant_species = levels(df_single$ant_species),
  forest_type = "primary",
  prey_shape  = mean(df_single$prey_shape, na.rm = TRUE))

# LRM-Predictions (rel_load3_l)
pred_rel_load_LRM <- predict(
  rel_load3_l,
  newdata  = newdat_rel_load_l,
  interval = "confidence")

newdat_rel_load_l$fit_log_rel_load <- pred_rel_load_LRM[, "fit"]
newdat_rel_load_l$lwr_log_rel_load <- pred_rel_load_LRM[, "lwr"]
newdat_rel_load_l$upr_log_rel_load <- pred_rel_load_LRM[, "upr"]

# Plot with regression line and confidence bands 
plot_rel_load_lengthLM <- ggplot(df_single, aes(x = ant_length,
                                                y = relativ_loading,
                                                colour = ant_species)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_line(data = newdat_rel_load_l,aes(x = ant_length,
                                         y = 10^fit_log_rel_load,
                                         colour = ant_species),
            linewidth = 1,
            inherit.aes = FALSE) +
  geom_ribbon( data = newdat_rel_load_l, aes( x = ant_length,
                                              ymin = 10^lwr_log_rel_load,
                                              ymax = 10^upr_log_rel_load,
                                              fill = ant_species),
               alpha = 0.4,
               colour = NA,
               inherit.aes = FALSE) +
  scale_x_log10(name = "Ant length (mm)") +
  scale_y_log10(name = "Relative loading (prey weight / ant weight)") +
  labs(colour = "Dorylus", fill = "Dorylus") +
  lrm_style +
  lrm_legend_theme

plot_rel_load_lengthLM

jpeg(file = "Relative Load vs. Ant Length LRM.jpg",
     width = 20, height = 18, units = "cm", res = 300)
plot_rel_load_lengthLM
dev.off()


# Residualplots for rel_load3_l
jpeg(file = "Residualplots Relative load LRM length.jpg",
     width = 20, height = 18, units = "cm", res = 300)
par(mfrow = c(2, 2))
plot(rel_load3_l)
dev.off()
#=============================================================================
                
##### 3.4 Loading in ALL Workers  #### 
#Load per ant vs ant weight + single vs multiple carriers
#=============================================================================#
               ######## A.  Model calculations ####

# It´s calculated for single and multiple carriers
df["ant_loading"] <- df$ant_weight / df$prey_weight

#relativ loading 
df["relativ_loading"] <- df$ant_loading / (1000*df$ant_weight)

#load per ant to get ant carring contribution in multiple carriers
df["load_per_ant"] <- df$ant_loading / df$total_ant_number

# categorize carrier state in each sample id 
df$carrier_type <- ifelse(df$total_ant_number > 1, "multiple", "single")
df$carrier_type <- factor(df$carrier_type)


#Load per ant in multiple vs. single carriers 
#herre we want to test if the diferent ant species behave differntly in their
#load per ant between carrier typs 

#Model 1 
carrier_model_1 <- lm(log10(load_per_ant) ~ log10(1000*ant_weight) + carrier_type 
                      + ant_species, data = df)
summary(carrier_model_1)
confint(carrier_model_1)
bptest(carrier_model_1)
coeftest(carrier_model_1, vcov = vcovHC(carrier_model_1, type = "HC3"))

#Model 2
carrier_model_2 <- lm(log10(load_per_ant) ~ log10(1000*ant_weight) + carrier_type 
                      * ant_species, data = df)
summary(carrier_model_2)
confint(carrier_model_2)
bptest(carrier_model_2)
coeftest(carrier_model_2, vcov = vcovHC(carrier_model_2, type = "HC3"))

#Model 3 
carrier_model_3 <- lm(log10(load_per_ant) ~ log10(1000*ant_weight) + carrier_type 
                    * ant_species + forest_type + prey_shape, data = df)

summary(carrier_model_3)
confint(carrier_model_3)
bptest(carrier_model_3)
coeftest(carrier_model_3, vcov = vcovHC(carrier_model_3, type = "HC3"))

#model 4 <- use this one
carrier_model_4 <- lm(log10(load_per_ant) ~ log10(1000*ant_weight) * carrier_type 
                           * ant_species + forest_type + prey_shape, data = df)

summary(carrier_model_4)
confint(carrier_model_4)
bptest(carrier_model_4)
coeftest(carrier_model_4, vcov = vcovHC(carrier_model_4, type = "HC3"))
anova(carrier_model_3, carrier_model_4)

jpeg(file = "Residualplots load of carrier typ.jpg",
     width = 20, height = 18, units = "cm", res = 300)
par(mfrow = c(2, 2))
plot(carrier_model_4)
dev.off()


AIC(carrier_model_1, carrier_model_2, 
    carrier_model_3, carrier_model_4)
#anova(carrier_model_1, carrier_model_2, carrier_model_3, carrier_model_4)
#anova not usefull since models not nested 


               ######## B.  LM Plot ####
#=============================================================================#
#setting up sequence 
seq_load_per <- seq(
  from = min(df$ant_weight, na.rm = TRUE),
  to   = max(df$ant_weight, na.rm = TRUE),
  length.out = 100
)

#Prediction grid (both species, primary forest, mean prey shape)
#Note: kick out species, forest, shpae. make "new model" just for weidht and 
#carrier typ and visualize this one. still talk about it in text 
newdat_load_per <- expand.grid(
  ant_weight     = seq_load_per,
  carrier_type = levels(df$carrier_type),
  ant_species  = levels(df$ant_species),
  forest_type  = "secondary",
  prey_shape   = mean(df$prey_shape, na.rm = TRUE)
  )

#Fixed-effect LM predictions
pred_load_per <- predict(carrier_model_3,
                                         newdata = newdat_load_per,
                                         interval = "confidence")

newdat_load_per$fit_log  <- pred_load_per[, "fit"]
newdat_load_per$lwr_log  <- pred_load_per[, "lwr"]
newdat_load_per$upr_log  <- pred_load_per[, "upr"]


load_per_mod <- ggplot(df, aes(x = 1000 * ant_weight,
                               y = load_per_ant)) +
  # Punkte: Farbe = carrier_type, Form = ant_species
  geom_point(aes(color = carrier_type,
                 shape = ant_species),
             alpha = 0.4,
             size  = 1) +
  
  # Modell-Linien: gleiche Farbzuordnung, Linientyp = ant_species
  geom_line(data = newdat_load_per,
            aes(x        = 1000 * ant_weight,
                y        = 10^fit_log,
                color    = carrier_type,
                linetype = ant_species,
                group    = interaction(carrier_type, ant_species)),
            linewidth = 0.8) +
  scale_x_log10(name = "Ant weight [mg]") +
  scale_y_log10(name = "Load per ant [mg]") 
  
  

load_per_mod

#==================================.===========================================

                        #### 4 Transport of Prey  #### 

#=============================================================================# 
##### 4.1 Frequency s of Carries #####

#creating a matrix with % of carriers
carrier_frequencies <- df_dis %>%
  group_by(ant_species, carrier_type) %>%
  summarise(N = n(), .groups = "drop") %>%
  group_by(ant_species) %>%
  mutate(Percent = round(100 * N / sum(N), 2))

carrier_frequencies

tab <- table(df_dis$ant_species, df_dis$carrier_type)
chisq.test(tab)

