#=============================================================================#
#         Statistical Analysis for the Dorylus prey Spectra Relations         #
#=============================================================================#
#### Packages and Version Controll####
## Code runs ons R.4.3.3 
## Some used packages are no longer supported since R.4.4.0 
## To date 18.11.25 no patch for R4.5. loading code will crash R

#Compendency for cooccur
install.packages("gmp")
library("gmp")
#load cooccour and EcoSimR directly from HATTPS since no longer in CRAN
install.packages(
  "https://cran.r-project.org/src/contrib/Archive/cooccur/cooccur_1.3.tar.gz",
  repos = NULL, type = "source") # cooccur
install.packages(
  "https://cran.r-project.org/src/contrib/Archive/EcoSimR/EcoSimR_0.1.0.tar.gz",
  repos = NULL, type = "source") 

install.packages("tidyverse")
install.packages("devtools")
install_github("fawda123/ggord")
install.packages("data.table")
install.packages("ggthemes")
install.packages("RColorBrewer")
install.packages("ggord")
install.packages("EcoSimR")
install.packages("vegan3d")
install.packages("rstatix")
install.packages("ggpubr")
install.packages("plot3D")
install.packages("ggtern")
install.packages("ggthemes")
install.packages("EcoSimR")
install.packages("mvabund")
install.packages("cooccur")
install.packages("ggord")
install.packages("viridis")
install.packages("ggforce")
install.packages("tidyquant")
install.packages("bipartite")
install.packages("Hmisc") # for mean_sdl
install.packages("ggbeeswarm")
install.packages("ggdist")

#### Library ####
library(gmp)                                #⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⢿⣿⣿⡿⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿
library(bipartite)                          #⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⢀⣾⣿⣿⣷⡀⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿
library(ggforce)                            #⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣇⠘⣿⣿⣿⣿⠃⣸⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿
library(tidyquant)                          #⣿⣿⣿⣿⣏⣉⠙⠛⠿⣿⣿⣿⣦⡾⠋⠙⢷⣴⣿⣿⣿⠿⠛⠋⣉⣹⣿⣿⣿⣿
library(EcoSimR)                            #⣿⣿⣿⣿⣿⣿⣿⣷⣦⡀⢻⣿⣟⠀⠀⠀⠀⣻⣿⡟⢀⣴⣾⣿⣿⣿⣿⣿⣿⣿
library(vegan)                              #⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⣄⠙⢿⣷⢤⡤⣾⡿⠋⣠⣾⣿⣿⣿⣿⣿⣿⣿⣿⣿
library(ggeffects)                          #⣿⣿⣿⣿⣿⣿⣿⣿⠉⠙⠛⠷⢼⡇⠀⠀⢸⡧⠾⠛⠋⠉⣿⣿⣿⣿⣿⣿⣿⣿
library(mvabund)                            #⣿⣿⣿⣿⣿⣿⣿⠇⢰⣿⣿⣶⡼⣷⠀⠀⣾⢧⣶⣾⣿⡆⠸⣿⣿⣿⣿⣿⣿⣿
library(ggplot2)                            #⣿⣿⣿⣿⣿⣿⡏⢠⣿⣿⠿⠋⣠⣿⡶⢶⣿⣄⠙⠿⣿⣿⡄⢹⣿⣿⣿⣿⣿⣿
library(effectsize)                         #⣿⣿⣿⣿⣿⡟⢠⣾⣿⠁⣴⣾⠋⠁⠀⠀⠈⠙⣷⣦⠈⣿⣷⡄⢻⣿⣿⣿⣿⣿
library(cooccur)                            #⣿⣿⣿⣿⠏⣠⣿⣿⣿⠀⣿⣧⠀⠀⠀⠀⠀⠀⣼⣿⠀⣿⣿⣿⣄⠹⣿⣿⣿⣿
library(tidyverse)                          #⣿⣿⣿⣷⣼⣿⣿⣿⡟⢀⣿⣿⡄⠀⠀⠀⠀⢠⣿⣿⡀⢻⣿⣿⣿⣧⣾⣿⣿⣿
require(graphics)                           #⣿⣿⣿⣿⣿⣿⣿⣿⠇⣸⣿⣿⣿⣄⠀⠀⣠⣿⣿⣿⣇⠸⣿⣿⣿⣿⣿⣿⣿⣿
library(devtools)                           #⣿⣿⣿⣿⣿⡿⠛⣁⣴⣿⣿⣿⣿⣿⣦⣴⣿⣿⣿⣿⣿⣦⣈⠛⢿⣿⣿⣿⣿⣿
library(ggord)                              #⣿⣿⣿⣿⣿⣷⣾⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⣾⣿⣿⣿⣿⣿
library(ggthemes)
library(RColorBrewer)
library(viridis)
library(ggpubr)
library(rstatix)
library(plot3D)
library(patchwork)
library(dplyr)
library(gridExtra)
library(SpadeR)
library(iNEXT)
library(Hmisc)
library(lme4)
library(ggbeeswarm)
library(ggdist)
library(MASS)
library(lmtest)
library(performance)
library(DHARMa)
library(lme4) #GLMM




#-----------------------------------------------------------------------------#
                    #### 1.Diet Composition by Ant Species ####
#-----------------------------------------------------------------------------#
# Data load in and preparation
readin <- read.csv("pooled prey.csv", header = TRUE) 
ants <- readin %>% remove_rownames %>% column_to_rownames(var = "raid_ID")

#genereate mulitvariant response with prey counts for manyglm
mvprey <- mvabund(ants[,3:27])

##### 1.1 Fitting GLMs for Multivariate Abundance Data ####
mod1 <- manyglm(mvprey ~ ants$ant_species * ants$forest_type,
                family = "negative.binomial",
                data = ants)

mvabund.aov <- anova.manyglm(mod1, p.uni = "adjusted",
                             resamp = "montecarlo", test = "LR")
#Results 
#First table shows global MGLMM results
#Furterh each toxon is includes individually 
mvabund.aov 

#Saving results in spreadsheets
#Note: in the univariate test, probably only Coleoptera significant for us
write.table(mvabund.aov$uni.p, "mod1 p values.csv")
write.table(mvabund.aov$uni.test, "mod1 test scores.csv")

##### 1.2 NMDS on diet composition by forest type and species #####

#Data load in and preparation
nmds1 <- read.csv("pooled prey.csv")
prey <- nmds1[,4:28] # this was done before
rownames(prey) <- nmds1$raid_ID

prey1.mds <- metaMDS(comm = prey, distance = "binomial", trymax=100, k=2, 
                     trace = FALSE, autotransform = T, 
                     expand = T, zerodist = "add")
summary(prey1.mds)
plot(prey1.mds$points)
identify(prey1.mds$points) 
prey1.mds$stress 
#NOTE:
#Lower the better, great <0.1, ideally <0.2. If >0.2 include additional dimension
#consider transforming matrix pre-MDS to reduce the influence of extreme values

# NMDS by ant species (colours) and forest type (shapes)
points <- prey1.mds$points
points <- as.data.frame(points)
mdsarmy3 <- ggplot(points, aes(x = MDS1, y = MDS2,
                               col = nmds1$ant_species,
                               shape = nmds1$forest_type)) + 
  geom_point(size = 3, alpha = 0.8) + 
  stat_ellipse(data = points, aes(MDS1, MDS2, group = c(nmds1$ant_species)),
               geom = "polygon", alpha=0.1, level = 0.80) +
  stat_ellipse(data = points, aes(MDS1, MDS2, group = c(nmds1$ant_species),
                                  col = c(nmds1$ant_species)), level = 0.8) +
  scale_color_manual(values = c("firebrick", "steelblue")) +
  scale_shape_manual(values = c("primary" = 16, "secondary" = 17)) +
  labs(x = "", y = "", fill = "Species") + theme_classic(base_size = 17) +  
  theme(legend.position = "right") +
  guides(color=guide_legend(title="Species")) + 
  guides(shape=guide_legend(title = "Forest type"))
mdsarmy3

# filesave
jpeg(file = "point cloud NMDS.jpg",
     width = 25, height = 18, units = "cm", res = 300)
mdsarmy3
dev.off()

##### 1.3 Pianka niche overlap #####
# This analysis is supplementary to the composition analysis
# Both measure whether the two species are consuming different things.

spdata <- read.csv("pianka dorylus.csv", header = TRUE) 
# Disregard the warning message!!!!

spmod <- niche_null_model(spdata, algo = "ra3", 
                          metric = "pianka", nReps = 10000)
summary(spmod) 
# Extremely high diet overlap levels. Observed Index = 0.95, represents almost 
# complete overlap. If both species were feeding at random, overlap = 0.26.
# Because the observed 0.95 value is so much higher than the simulated 0.26,
# very clear sign that the niche is almost totally overlapping.

plot(spmod, type = "hist") 
# the plotting function kinda sucks. Don't include in papers.

##### 1.4 Bipartite plots #####
# Convert Diet frequency to proportional data
bpno <- read.csv("pianka dorylus prop.csv") 
summary(bpno)
rownames(bpno) <- bpno[,1]
rownames(bpno)
colnames(bpno)

plotweb(bpno[,2:26], high.spacing = 0.005, low.spacing = 0.1,
        arrow = "up.centre", text.rot = 90, bor.col.high = "black", 
        method = "normal", 
        col.interaction = viridis(25),
        bor.col.interaction = viridis(25),
        labsize = 1,  plot.axes = F, high.y = 1.5, low.y = 0.6)

jpeg(file = "dorylus diet bipartite.jpg", width = 25, height = 20, 
     units="cm", res=300)
plotweb(bpno[,2:26], high.spacing = 0.005, low.spacing = 0.1,
        arrow = "up.centre", text.rot = 90, bor.col.high = "black", 
        method = "normal", 
        col.interaction = viridis(25),
        bor.col.interaction = viridis(25),
        labsize = 1,  plot.axes = F, high.y = 1.5, low.y = 0.6)
dev.off()

#-----------------------------------------------------------------------------#
            #### 2.Dismembered - Preparation for transport####
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
#analysis where made with ant species and forest typ as factor

# GLM on Multivariate Abundance Data (mvabund)
lifestage <- read.csv("prey spectra lifestage.csv", header = T)

# Count number of prey in each life stage per ant species
life_counts <- lifestage %>%
  group_by(raid_ID, ant_species, forest_type, life_stage) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = life_stage, values_from = count, values_fill = 0)
life_counts

mvprey <- mvabund(life_counts[,4:8])

#NOTE: Still need to make it a bit more fancy and look into whats happening here
mod1 <- manyglm(mvprey ~ life_counts$ant_species * life_counts$forest_type,
                family = "negative.binomial",
                data = life_counts)
plot.manyglm(mod1)
plot(mod1)

mvabund.aov <- anova.manyglm(mod1, p.uni = "adjusted", 
                             resamp = "montecarlo", test = "LR")
mvabund.aov

write.table(mvabund.aov$uni.p, "mod1 p values.csv")
write.table(mvabund.aov$uni.test, "mod1 test scores.csv")


