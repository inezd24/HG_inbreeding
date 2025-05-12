# Script developed by Inez Derkx, summer 2024
# Script aims to examine the relationship between ROH-based inbreeding and subsistence style

##################################### IMPORT LIBRARIES #####################################

library(ggplot2)
library(reshape2)
library(readr)
library(dplyr)
library(brms)
library(cmdstanr)

##################################### IMPORT DATA #####################################

# Import necessary files
populations <- c("BaYaka", "Biaka", "Agta","Cameroon_Mbo", "Controls", "Dailekh", "Datog", "Hadza", "Raute","Sandawe")
for (pop in populations) {
  assign(paste(pop, ".roh", sep = ""), read.csv(paste("/location/of/files/", pop, ".hom", sep = ""), sep = ""))
}

# Rename these populations
Controls.roh$FID <- 'Palanan'
Cameroon_Mbo.roh$FID <- 'Mbo'

# Combine files
roh_farm_forage <- rbind(BaYaka.roh, Biaka.roh, Agta.roh, Cameroon_Mbo.roh, Controls.roh, Dailekh.roh, Datog.roh, Hadza.roh, Raute.roh, Sandawe.roh)

# Categorize by subsistence
roh_farm_forage$Subsistence <- ifelse(roh_farm_forage$FID == 'Agta' | roh_farm_forage$FID == 'Raute' | roh_farm_forage$FID == 'BaYaka' |
                                      roh_farm_forage$FID == 'Biaka' | roh_farm_forage$FID == 'Mbuti' | roh_farm_forage$FID == 'Hadza', 'Hunter-Gatherer', 'Other subsistence')

# Order population name according to region
roh_farm_forage$FID <- factor(roh_farm_forage$FID, levels = c("Agta", "Palanan", "BaYaka", "Biaka", "Mbo",  "Hadza", "Datog", "Sandawe", "Raute", "Dailekh"))

# Add region (continent)
roh_farm_forage$region <- ifelse(roh_farm_forage$FID == 'Agta' | roh_farm_forage$FID == 'Palanan', 'Southeast Asia',
                                 ifelse(roh_farm_forage$FID == 'Raute' | roh_farm_forage$FID == 'Dailekh', 'South Asia',
                                        ifelse(roh_farm_forage$FID == 'Biaka' | roh_farm_forage$FID == 'BaYaka' | roh_farm_forage$FID == 'Mbo', 'West Africa', 'East Africa')))

# Order region according to proximity
roh_farm_forage$region <- factor(roh_farm_forage$region, levels = 
                                   c("East Africa", "West Africa", "South Asia", "Southeast Asia"))

# Consider only segments longer than 1.5 MB
roh_farm_forage <- roh_farm_forage[(roh_farm_forage$KB > 1500),]
roh_farm_forage$MB <- roh_farm_forage$KB/1000

# Compute SROH, NROH, and FROH
roh_farm_forage <- roh_farm_forage %>% group_by(Subsistence, region, FID, IID) %>% summarise(SROH = sum(KB), NROH = n(), FROH = sum(KB)/2881033)

# Remove Datog and Biaka > redundant in dataset due to population size / other constraints
roh_farm_forage <-roh_farm_forage[!(roh_farm_forage$FID == 'Datog'),]
roh_farm_forage <-roh_farm_forage[!(roh_farm_forage$FID == 'Biaka'),]

##################################### DESCRIPTIVE VISUALISATION OF DATA #####################################

# Examine relationship between inbreeding (FROH), subsistence style, and region
# MANUSCRIPT FIGURE 1
png("/location/of/plot/FROH.png", width = 1600, height = 800)
ggplot(roh_farm_forage, aes(x = FID, y = FROH, fill = Subsistence)) +
  geom_violin(position = position_dodge(1), trim = F, alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white") +
  scale_fill_manual(values = c("#ffa635", "#839ffe")) +
  facet_wrap(~region, nrow=1, ncol = 4, scales = "free_x") +  # Facet by region with independent x-scales
  theme_classic() +
  theme(text = element_text(size = 30),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  labs(y = "FROH", x = "Population", fill = 'Subsistence style')
dev.off()



##################################### MODELLING OF DATA #####################################

## Base model: just FROH and Subsistence
# Perform bayesian model
froh_subsistence <- brm(FROH ~ Subsistence,
                        data = roh_farm_forage,
                        family = Beta(),
                        prior = c(
                          prior(normal(0, 1), class = "b")  # default weakly informative
                        ),
                        backend = "cmdstanr")

# Examine the posterior distribution
summary(froh_subsistence)
plot(froh_subsistence, variable = c("b_SubsistenceOthersubsistence"))
cond_effects1 <- conditional_effects(froh_subsistence, effects = "Subsistence")
ce1 <- conditional_effects(froh_subsistence)
cedata1 <- ce1$Subsistence # for group means and CIs
cond_data1 <- cond_effects1[[1]]  

# Plot the conditional effects of subsistence on FROH
png("/Users/inezd/Documents/Science/PhD/Chapter_3/Plots/FROH_subsistence.png", width = 800, height = 800)
ggplot(cond_data1, aes(x = Subsistence, y = estimate__, ymin = lower__, ymax = upper__, group = Subsistence)) +
  geom_pointrange(linewidth = 1) +  # Plot the point estimate with confidence intervals
  geom_point(aes(color = Subsistence),size = 5, alpha = 0.6) +  # Customize point appearance
  theme_classic() +
  scale_color_manual(values=c("#ffa635", "#839ffe")) +
  theme(text = element_text(size = 25),
        legend.position = "none") +
  labs(y = "Estimated FROH", x = "Subsistence")
dev.off()


### MODEL 1: Subsistence and FROH, controlled by region

# Check the model prerequisites
hist(roh_farm_forage$FROH) #NOT normally distributed and since 0 < FROH < 1, we use 'Beta' as family

# Perform bayesian model
froh_subsistence_region <- brm(FROH ~ Subsistence + (1|region),
                               family = Beta(),
                               data = roh_farm_forage,
                               prior = c(
                                 prior(normal(0, 2), class = "b"),  # Regularizing prior on fixed effects
                                 prior(student_t(3, 0, 2.5), class = "Intercept"),  # Weakly informative prior on intercept
                                 prior(cauchy(0, 1), class = "sd")  # Prior on group-level standard deviation
                               ),
                               control = list(adapt_delta = 0.95), # Useful for models with complex priors
                               backend = "cmdstanr")  

# Examine the posterior distribution
summary(froh_subsistence_region)
plot(froh_subsistence_region, variable = c("b_SubsistenceOthersubsistence"))
cond_effects <- conditional_effects(froh_subsistence_region, effects = "Subsistence")
ce <- conditional_effects(froh_subsistence_region)
cedata <- ce$Subsistence # for group means and CIs
cedata
cond_data <- cond_effects[[1]]  

# Plot the conditional effects of subsistence on FROH
# Supplementary Figure 1a
png("/Users/inezd/Documents/Science/PhD/Chapter_3/Plots/FROH_subsistence_region.png", width = 800, height = 800)
ggplot(cond_data, aes(x = Subsistence, y = estimate__, ymin = lower__, ymax = upper__, group = Subsistence)) +
  geom_pointrange(linewidth = 1) +  # Plot the point estimate with confidence intervals
  geom_point(aes(color = Subsistence),size = 5, alpha = 0.6) +  # Customize point appearance
  theme_classic() +
  scale_color_manual(values=c("#ffa635", "#839ffe")) +
  theme(text = element_text(size = 25),
        legend.position = "none") +
  labs(y = "Estimated FROH", x = "Subsistence")
dev.off()



### MODEL 2: FROH and interaction Subsistence and Region

# Perform bayesian model
froh_int_subsistence_region <- brm(
  FROH ~ Subsistence * region,
  data = roh_farm_forage,
  family = Beta(),
  prior = c(
    prior(normal(0, 1), class = "b"),              # Priors for fixed effects
    prior(normal(0, 1), class = "Intercept"),      # Prior for intercept
    prior(exponential(1), class = "phi")           # Prior for precision (phi) in Beta
  ),
  backend = "cmdstanr"
)

# Examine the posterior distribution
summary(froh_int_subsistence_region)
plot(froh_int_subsistence_region, variable = c("b_SubsistenceOthersubsistence"))
ce_int <- conditional_effects(froh_int_subsistence_region)
cedata <- ce_int$`Subsistence:region` # for group means and CIs
cedata
cond_effects_int <- conditional_effects(froh_int_subsistence_region, effects = "Subsistence:region")
cond_data_int <- cond_effects_int[[1]]  

# Plot the conditional effects of subsistence on FROH
# Supplementary Figure 1b
png("/Users/inezd/Documents/Science/PhD/Chapter_3/Plots/FROH_int_subsistence_region.png", width = 800, height = 800)
ggplot(cond_data_int, aes(x = region, y = estimate__, ymin = lower__, ymax = upper__, color = Subsistence, group = Subsistence)) +
  geom_pointrange(position = position_dodge(width = 0.4), linewidth = 1) +
  geom_point(position = position_dodge(width = 0.4), size = 5, alpha = 0.6) +
  theme_classic() +
  scale_color_manual(values = c("#ffa635", "#839ffe")) +
  theme(text = element_text(size = 25)) +
  labs(y = "Estimated FROH", x = "Region", color = "Subsistence")
dev.off()


## end of script ##
