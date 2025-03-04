# Script developed by Inez Derkx, summer 2024
# Script aims to examine the relationship between ROH-based inbreeding and subsistence style

##################################### IMPORT LIBRARIES #####################################

library(ggplot2)
library(reshape2)
library(readr)
library(dplyr)
library(brms)

##################################### IMPORT DATA #####################################

# Import necessary files
roh_compare <- read.csv("/Users/inezd/PHD/Genetic_Data/homozyg/roh_compare.csv", sep = ',')
roh_compare$X <- NULL
populations <- c("BaYaka", "Biaka", "Agta","Cameroon_Mbo", "Controls", "Dailekh", "Datog", "Hadza", "Raute","Sandawe")
for (pop in populations) {
  assign(paste(pop, ".roh", sep = ""), read.csv(paste("/Users/inezd/PHD/Genetic_Data/homozyg/RUNSoh/", pop, ".hom", sep = ""), sep = ""))
}
Controls.roh$FID <- 'Palanan'
Cameroon_Mbo.roh$FID <- 'Mbo'

# Transform files
roh_farm_forage <- rbind(BaYaka.roh, Biaka.roh, Agta.roh, Cameroon_Mbo.roh, Controls.roh, Dailekh.roh, Datog.roh, Hadza.roh, Raute.roh, Sandawe.roh)
roh_farm_forage$Subsistence <- ifelse(roh_farm_forage$FID == 'Agta' | roh_farm_forage$FID == 'Raute' | roh_farm_forage$FID == 'BaYaka' |
                                      roh_farm_forage$FID == 'Biaka' | roh_farm_forage$FID == 'Mbuti' | roh_farm_forage$FID == 'Hadza', 'Hunter-Gatherer', 'Other subsistence')
roh_farm_forage$FID <- factor(roh_farm_forage$FID, levels = c("Agta", "Palanan", "BaYaka", "Biaka", "Mbo",  "Hadza", "Datog", "Sandawe", "Raute", "Dailekh"))
roh_farm_forage$region <- ifelse(roh_farm_forage$FID == 'Agta' | roh_farm_forage$FID == 'Palanan', 'Southeast Asia',
                                 ifelse(roh_farm_forage$FID == 'Raute' | roh_farm_forage$FID == 'Dailekh', 'South Asia',
                                        ifelse(roh_farm_forage$FID == 'Biaka' | roh_farm_forage$FID == 'BaYaka' | roh_farm_forage$FID == 'Mbo', 'West Africa', 'East Africa')))
roh_farm_forage$region <- factor(roh_farm_forage$region, levels = 
                                   c("East Africa", "West Africa", "South Asia", "Southeast Asia"))

# Consider only segments longer than 1.5 MB
roh_farm_forage <- roh_farm_forage[(roh_farm_forage$KB > 1500),]
roh_farm_forage$MB <- roh_farm_forage$KB/1000
roh_farm_forage <- roh_farm_forage %>% group_by(Subsistence, region, FID, IID) %>% summarise(SROH = sum(KB), NROH = n(), FROH = sum(KB)/2881033)

# Remove Datog and Biaka > they only have a three sample dataset
roh_farm_forage <-roh_farm_forage[!(roh_farm_forage$FID == 'Datog'),]
roh_farm_forage <-roh_farm_forage[!(roh_farm_forage$FID == 'Biaka'),]


##################################### DESCRIPTIVE VISUALISATION OF DATA #####################################


# Examine relationship between inbreeding (FROH), subsistence style, and region
png("/Users/inezd/PHD/Chapter_3/Plots/FROH_by_subsistence.png", width = 1600, height = 800)
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

# For supplementary info, do the same for the NROH and SROH
png("/Users/inezd/PHD/Chapter_3/Plots/SROH_by_subsistence.png", width = 800, height = 800)
ggplot(roh_farm_forage, aes(x = FID, y = SROH, fill = Subsistence)) +
  geom_violin(position = position_dodge(1), trim = F, alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white") +
  scale_fill_manual(values = c("#ffa635", "#839ffe")) +
  facet_wrap(~region, nrow=1, ncol = 4, scales = "free_x") +  # Facet by region with independent x-scales
  theme_classic() +
  theme(text = element_text(size = 30),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  labs(y = "SROH", x = "Population", fill = 'Subsistence style')
dev.off()

# Plot NROH
png("/Users/inezd/PHD/Chapter_3/Plots/NROH_by_subsistence.png", width = 800, height = 800)
ggplot(roh_farm_forage, aes(x = FID, y = NROH, fill = Subsistence)) +
  geom_violin(position = position_dodge(1), trim = F, alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white") +
  scale_fill_manual(values = c("#ffa635", "#839ffe")) +
  facet_wrap(~region, nrow=1, ncol = 4, scales = "free_x") +  # Facet by region with independent x-scales
  theme_classic() +
  theme(text = element_text(size = 30),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  labs(y = "NROH", x = "Population", fill = 'Subsistence style')
dev.off()


##################################### MODELLING OF DATA #####################################

### MODEL 1: Subsistence and FROH, controlled by region

# Check the model prerequisites
hist(roh_farm_forage$FROH) #NOT normally distributed and since 0 < FROH < 1, we use 'Beta' as family

# Perform bayesian model
froh_subsistence <- brm(FROH ~ Subsistence + (1|region),
            family = Beta(),
            data = roh_farm_forage,
            prior = c(
              prior(normal(0, 2), class = "b"),  # Regularizing prior on fixed effects
              prior(student_t(3, 0, 2.5), class = "Intercept"),  # Weakly informative prior on intercept
              prior(cauchy(0, 1), class = "sd")  # Prior on group-level standard deviation
            ),
            control = list(adapt_delta = 0.95)  # Useful for models with complex priors
)

# Examine the posterior distribution
summary(froh_subsistence)
plot(froh_subsistence, variable = c("b_SubsistenceOthersubsistence"))
cond_effects <- conditional_effects(froh_subsistence, effects = "Subsistence")
cond_data <- cond_effects[[1]]  

# Plot the conditional effects of subsistence on FROH
png("/Users/inezd/PHD/Chapter_3/Plots/FROH_subsistence.png", width = 800, height = 800)
ggplot(cond_data, aes(x = Subsistence, y = estimate__, ymin = lower__, ymax = upper__, group = Subsistence)) +
  geom_pointrange(linewidth = 1) +  # Plot the point estimate with confidence intervals
  geom_point(aes(color = Subsistence),size = 5, alpha = 0.6) +  # Customize point appearance
  theme_classic() +
  scale_color_manual(values=c("#ffa635", "#839ffe")) +
  theme(text = element_text(size = 25),
        legend.position = "none") +
  labs(y = "Estimated FROH", x = "Subsistence")
dev.off()



### MODEL 2: Subsistence, population size, and FROH, controlled by region

# Add meta population size based on published records
roh_farm_forage$popsize <- ifelse(roh_farm_forage$FID == "Agta", 10000, 
                               ifelse(roh_farm_forage$FID == 'Raute', 150, 
                                      ifelse(roh_farm_forage$FID == 'Hadza', 1000, 
                                                    ifelse(roh_farm_forage$FID == 'BaYaka', 190000, 
                                                                  ifelse(roh_farm_forage$FID == 'Sandawe', 30000, NA))))) # Sandawe: 30000 (Lachance, 2012)
roh_foragers <- roh_farm_forage[(roh_farm_forage$Subsistence == 'Hunter-Gatherer'),]

# Models will not converge if we don't have population size for all populations. 


## end of script ##
