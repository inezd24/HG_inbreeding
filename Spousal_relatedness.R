# Script developed by Inez Derkx, summer 2024
# Script aims to examine list of spouses and their relatedness alongside relatedness between all those 22 individuals. 


##################################### IMPORT LIBRARIES #####################################

library(ggplot2)
library(reshape2)
library(readr)
library(dplyr)

##################################### IMPORT DATA #####################################

raute_census_combined <- readxl::read_xlsx("/Users/inezd/Hunter Gatherer Resilience/Input/raute_census_january_2022.xlsx",
                                           sheet = "Census combined")
raute_only_IBD <- read_table("/Users/inezd/PHD/Genetic_Data/IBD_analyses/all.raute.IBD.Merged", 
                             col_names = c("sample_1", "haplo_1", "sample_2", "haplo_2", "chr", "start", "end","LOD", "length_cM"))
raute_genealogies_complete <- readxl::read_xlsx("/Users/inezd/Hunter Gatherer Resilience/Input/raute_census_january_2022.xlsx",
                                                sheet = "genealogies - correct")
raute_20_kin <- read.delim("~/PHD/Genetic_Data/IBIS/raute_maf_20.coef", 
                           col.names = c("sample_1", "sample_2", "kinship_20", "IBD2_20", "segment_count", "degree_of_relatedness"))
raute_7_kin <- read.delim("~/PHD/Genetic_Data/IBIS/raute_maf_7.coef",
                          col.names = c("sample_1", "sample_2", "kinship_7", "IBD2_7", "segment_count", "degree_of_relatedness"))


# Add some new colours to the colour palette
own_colours = list(
  spring = c("#fd7f6f", "#7eb0d5", "#b2e061", "#bd7ebe", "#ffb55a", "#ffee65", "#beb9db", "#fdcce5", "#8bd3c7", 
             "#5aa3ff", "#f08080", "#b6fad3"),
  retro = c("#f46a9b", "#ef9b20", "#edbf33", "#ede15b", "#bdcf32", "#87bc45", "#27aeef", "#b33dc6"),
  dutch= c("#e60049", "#0bb4ff", "#50e991", "#e6d800", "#9b19f5", "#ffa300", "#dc0ab4", "#b3d4ff", "#00bfa0"),
  rivers = c("#b30000", "#7c1158", "#4421af", "#1a53ff", "#0d88e6", "#00b7c7", "#5ad45a", "#8be04e", "#ebdc78"),
  pastels = c("#54bebe", "#76c8c8", "#98d1d1", "#badbdb", "#dedad2", "#e4bcad", "#df979e", "#d7658b", "#c80064"),
  colors_paper = c("#A3CBFA", "#78A2D6", "#5080B0", "#034984", "#4D4D4D", "#B0B0B0",
                   "#F4E98C", "#F0B065", "#ED9434", "#C8751A", "#5B290C", "#000000",
                   "#F5DCE8", "#E099AA", "#CC6677", "#B0475D", "#764786", "#CD86E6",
                   "#EF048B", "#A0017C"),
  colourblind_friendly = c("#ffa635", "#b77fff", "#ff9dcc", "#839ffe"))


## SPOUSAL INFO

# Make list of spouses
spouses <- raute_genealogies_complete[!is.na(raute_genealogies_complete$partner_id) & (raute_genealogies_complete$partner_alive == 'yes') & (raute_genealogies_complete$alive == 'yes'),c(1,2,3,14,15)]
names(spouses)[c(1,2,5)] <- c("id_2014_1", 'id', 'old_id')
spouses <- merge(spouses, raute_genealogies_complete[,1:2], by = 'old_id', all.x=T)
names(spouses)[6] <- c("id_2014_2")
spouses <- spouses[c(4,5,2,6)]
spouses$pair <- paste(pmin(spouses$id_2014_1, spouses$id_2014_2), pmax(spouses$id_2014_2, spouses$id_2014_1), sep = '-')
spouses$id1 <- pmin(spouses$id_2014_1, spouses$id_2014_2)
spouses$id2 <- pmax(spouses$id_2014_1, spouses$id_2014_2)
spouses <- spouses[!duplicated(spouses$pair),6:7]

# Add oragene ids to spouses list
names(id_gen)[1:2] <- c('ind1', 'id1')
spouses <- merge(spouses, id_gen[,1:2], by = 'id1')
names(id_gen)[1:2] <- c('ind2', 'id2')
spouses <- merge(spouses, id_gen[,1:2], by = 'id2')
names(id_gen)[1:2] <- c("ind", "id_2014") # we have 11 couples included with the analysis 

# Merge the file with new ids >> make sure that you check the raute_fam file 
# column 8 should be new_IID, you can re do the file in fam_files script. 
names(raute_fam)[1] <- 'ind1'
spouses <- merge(spouses, raute_fam[,c(1,8)], by = 'ind1', all.x=T)
names(raute_fam)[1] <- 'ind2'
spouses <- merge(spouses, raute_fam[,c(1,8)], by = 'ind2', all.x=T)
names(raute_fam)[1] <- 'V1'
names(spouses)[5:6] <- c("sample_1", "sample_2")
spouses <- spouses[,c(2,1,4,3,5:6)]
spouses$pair <- paste(pmin(spouses$sample_1, spouses$sample_2), pmax(spouses$sample_1, spouses$sample_2), sep = '-')


## ADD IBD SHARING

# Summarize per 2 individuals
raute_only_IBD$pair <- paste(pmin(raute_only_IBD$sample_1, raute_only_IBD$sample_2), pmax(raute_only_IBD$sample_1, raute_only_IBD$sample_2), sep = '-')
raute_IBD_summary <- raute_only_IBD %>% group_by(pair) %>% summarise(n = n(), 
                                                                     sum = sum(length_cM), 
                                                                     mean = mean(length_cM))

# Summary stats
mean(raute_IBD_summary$sum) # Mean sum: 2207.67
mean(raute_IBD_summary$mean) # Mean mean: 14.91
mean(raute_IBD_summary$n) # Mean number: 169.17


# Merge this with the spousal info
spouses <- merge(spouses, raute_IBD_summary, by = c("pair"), all.x=T)

# Summary stats
mean(spouses$sum) # Mean sum: 2078.2
mean(spouses$mean) # Mean mean: 11.5
mean(spouses$n) # Mean number: 182.64


## ADD KINSHIP COEFFICIENTS 

# Add kinship coefficient: 20 cM
raute_20_kin$pair <- paste(pmin(raute_20_kin$sample_1, raute_20_kin$sample_2), pmax(raute_20_kin$sample_1, raute_20_kin$sample_2), sep = '-')
spouses <- merge(spouses, raute_20_kin[,c(7,3)], by = c("pair"), all.x=T)

# Add kinship coefficient: 7 cM
raute_7_kin$pair <- paste(pmin(raute_7_kin$sample_1, raute_7_kin$sample_2), pmax(raute_7_kin$sample_1, raute_7_kin$sample_2), sep = '-')
spouses <- merge(spouses, raute_7_kin[,c(7,3)], by = c("pair"), all.x=T)
spouses$type <- 'spouse'

# What are the means?
mean(spouses$kinship_20) # 0.07
mean(spouses$kinship_7) #0.12

# We have to check how these different variables correlate:
cor.test(spouses$sum, spouses$kinship_20) #0.0216, 0.68
cor.test(spouses$sum, spouses$kinship_7) #0.0009, 0.85 > so there is most overlap between these two
cor.test(spouses$n, spouses$kinship_20) #0.0264, -0.66
cor.test(spouses$n, spouses$kinship_7) #0.0294, -0.65

# Expected inbreeding offspring
expected_inbreeding_offspring <- function(father, mother, kinship) {
  return(0.5 * (father + mother) + 0.5 * kinship)
}


# What are the expected offspring inbreeding coefficients for these pairs?
names(roh_compare)[3] <- 'sample_1'
spouses <- merge(spouses, roh_compare[,c(3,6)], by = 'sample_1', all.x=T)
names(roh_compare)[3] <- 'sample_2'
spouses <- merge(spouses, roh_compare[,c(3,6)], by = 'sample_2', all.x=T)
names(roh_compare)[3] <- 'IID'
spouses$inbreeding_coef <- 0.5*(spouses$FROH.x + spouses$FROH.y) + (0.5*spouses$kinship_7)
mean(spouses$inbreeding_coef)



## EXAMINE AVERAGE AMONGST SPOUSES



# Now we want to examine: what is average the relatedness amongst these individuals?

# Create list of individuals who are spouses
spouse_ind <- reshape2::melt(spouses[,1:3], id.vars = 'pair')
spouse_ind <- data.frame(sample_1 = unique(spouse_ind$value))
spousal_relatedness <- data.frame(t(combn(spouse_ind$sample_1, 2)))
spousal_relatedness$pair <- paste(pmin(spousal_relatedness$X1, spousal_relatedness$X2), 
                                  pmax(spousal_relatedness$X1, spousal_relatedness$X2), 
                                  sep = '-')

# Add IBD list and kinship
spousal_relatedness <- merge(spousal_relatedness, raute_IBD_summary, by = 'pair', all.x=T)
spousal_relatedness <- merge(spousal_relatedness, raute_7_kin[,c(7,3)], by = c("pair"), all.x=T)
spousal_relatedness <- merge(spousal_relatedness, raute_20_kin[,c(7,3)], by = c("pair"), all.x=T)
spousal_relatedness <- merge(spousal_relatedness, spouses[,c(3,13)], by = 'pair', all.x=T)
spousal_relatedness$type[is.na(spousal_relatedness$type)] <- 'non-spouse'
table(spousal_relatedness$type) # 11 spouses, 220 non-spouses

# Summary stats
mean(spousal_relatedness$sum) # 2184.90
mean(spousal_relatedness$kinship_20) # 0.09
mean(spousal_relatedness$kinship_7)  # 0.14

observed_mean <- 0.124274  # Replace with your observed mean
t.test(spousal_relatedness$kinship_7, spouses$kinship_7) #0.006521

#Now: divided by spouse versus non-spouse
spousal_relatedness %>% group_by(type) %>% summarise(mean_IBD = mean(sum),
                                                     n_IBD = mean(n),
                                                 mean_7cM = mean(kinship_7), 
                                                 mean_20cM = mean(kinship_20)) 
# Spousal individuals have lower shared IBD and kinship coefficients 

# Actually check
anova_spouses <- aov(sum ~ type, data = spousal_relatedness)
summary(anova_spouses)
anova_spouses <- aov(kinship_7 ~ type, data = spousal_relatedness)
summary(anova_spouses)


# First examine just the distribution of all three values
png("/Users/inezd/PHD/Chapter_3/Plots/spousal_IBD_avg.png", height = 800, width = 800)
spousal_IBD <- ggplot(spousal_relatedness, aes(x = sum)) +
  geom_histogram(color = "#ebefff", fill="#839ffe", alpha = 0.6, bins = 50) +
  geom_segment(aes(x = 2207, xend = 2207, y = 0, yend = Inf), color = "#ffa635", linetype = 'dashed', linewidth = 1) + # orange // whole population 
  geom_segment(aes(x = 2078, xend = 2078, y = 0, yend = Inf), color = "#b77fff", linetype = 'dashed', linewidth = 1) + # amongst spouses // purple
  geom_segment(aes(x = mean(sum), xend = mean(sum), y = 0, yend = Inf), color =  "#ff9dcc", linetype = 'dashed', linewidth = 1) + #of this sample // pink
  theme_classic() +
  theme(text = element_text(size = 30)) +
  labs(title = "", x = 'Sum of IBD shared segments', y = 'frequency')
spousal_IBD
dev.off()
png("/Users/inezd/PHD/Chapter_3/Plots/spousal_kin7_avg.png", height = 800, width = 800)
spousal_kin7 <- ggplot(spousal_relatedness, aes(x = kinship_7)) +
  geom_histogram(color = "#ebefff", fill="#839ffe", alpha = 0.6, bins = 50) +
  geom_segment(aes(x = 0.14, xend = 0.14, y = 0, yend = Inf), color = "#ffa635", linetype = 'dashed', linewidth = 1) + # orange // whole population 
  geom_segment(aes(x = 0.12, xend = 0.12, y = 0, yend = Inf), color = "#b77fff", linetype = 'dashed', linewidth = 1) + # amongst spouses // purple
  geom_segment(aes(x = mean(kinship_7), xend = mean(kinship_7), y = 0, yend = Inf), color =  "#ff9dcc", linetype = 'dashed', linewidth = 1) + #of this sample // pink
  theme_classic() +
  theme(text = element_text(size = 30)) +
  scale_x_continuous(limits = c(0.05, 0.30), breaks = seq(0.05, 0.30, by = 0.05)) +
  labs(title = "", x = 'Kinship coefficient (7 cM)', y = 'frequency')
spousal_kin7
dev.off()
png("/Users/inezd/PHD/Chapter_3/Plots/spousal_kin20_avg.png", height = 800, width = 800)
spousal_kin20 <- ggplot(spousal_relatedness, aes(x = kinship_20)) +
  geom_histogram(color = "#ebefff", fill="#839ffe", alpha = 0.6, bins = 50) +
  geom_segment(aes(x = 0.09, xend = 0.09, y = 0, yend = Inf), color = "#ffa635", linetype = 'dashed', linewidth = 1) + # orange // whole population 
  geom_segment(aes(x = 0.07, xend = 0.07, y = 0, yend = Inf), color = "#b77fff", linetype = 'dashed', linewidth = 1) + # amongst spouses // purple
  geom_segment(aes(x = mean(kinship_20), xend = mean(kinship_20), y = 0, yend = Inf), color =  "#ff9dcc", linetype = 'dashed', linewidth = 1) + #of this sample // pink
  theme_classic() +
  theme(text = element_text(size = 30)) +
  scale_x_continuous(limits = c(0, 0.30), breaks = seq(0, 0.30, by = 0.05)) +
  labs(title = "", x = 'Kinship coefficient (20 cM)', y = 'frequency')
spousal_kin20
dev.off()



## DOES THIS MAKE ANY SENSE? > NOT REALLY



## whole sample

# Simulate sampling 1000 couples and calculate their IBD
set.seed(123) # For reproducibility
simulated_means_IBD <- replicate(1000, {
  sampled_pair <- spousal_relatedness[sample(nrow(spousal_relatedness), 11), ]
  return(mean(sampled_pair[["sum"]]))
})

# Simulate sampling 1000 couples and calculate their relatedness
set.seed(123) # For reproducibility
simulated_means_7cM <- replicate(1000, {
  sampled_pair <- spousal_relatedness[sample(nrow(spousal_relatedness), 11), ]
  return(mean(sampled_pair[["kinship_7"]]))
})

# Simulate sampling 1000 couples and calculate their relatedness
set.seed(123) # For reproducibility
simulated_means_20cM <- replicate(1000, {
  sampled_pair <- spousal_relatedness[sample(nrow(spousal_relatedness), 11), ]
  return(mean(sampled_pair[["kinship_20"]]))
})

#Compare these to the actual values of the 


# Plot distribution of simulated relatedness
simulated_relatednessIBD <- data.frame(sim = simulated_means_IBD)
mean(simulated_relatednessIBD$sim) #2183.916
mean(spouses$sum) #2078.23
png("/Users/inezd/PHD/Chapter_3/Plots/simulated_relatednessIBD.png", height = 800, width = 800)
ggplot(simulated_relatednessIBD, aes(x=sim)) +
  geom_histogram(color = "#ebefff", fill="#839ffe", alpha = 0.6, bins = 50) +
  geom_segment(aes(x = mean(spouses$sum), 
               xend = mean(spouses$sum), 
               y = 0, 
               yend = Inf), 
               color = "#b77fff", 
               linetype = 'dashed', 
               linewidth = 1) +  
  #ggtitle("Simulated kinship whole sample") +
  labs(x = "IBD sum length", y = 'Frequency') +
  theme_classic() +
  theme(text = element_text(size = 30))
dev.off()

# Plot distribution of simulated relatedness
simulated_relatedness7 <- data.frame(sim = simulated_means_7cM)
mean(simulated_relatedness7$sim) #0.1377262
mean(spouses$kinship_7) #0.124274
png("/Users/inezd/PHD/Chapter_3/Plots/simulated_relatedness7.png", height = 800, width = 800)
ggplot(simulated_relatedness7, aes(x=sim)) +
  geom_histogram(color = "#ebefff", fill="#839ffe", alpha = 0.6, bins = 50) +
  geom_segment(aes(x = mean(spouses$kinship_7), 
                   xend = mean(spouses$kinship_7), 
                   y = 0, 
                   yend = Inf), 
               color = "#b77fff", 
               linetype = 'dashed', 
               linewidth = 1) + 
  #ggtitle("Simulated kinship whole sample") +
  labs(x = "Relatedness (7cM)", y = 'Frequency') +
  theme_classic() +
  theme(text = element_text(size = 20))
dev.off()

# Plot distribution of simulated relatedness
simulated_relatedness20 <- data.frame(sim = simulated_means_20cM)
mean(simulated_relatedness20$sim) #0.08921228
mean(spouses$kinship_20) #0.07237264
png("/Users/inezd/PHD/Chapter_3/Plots/simulated_relatedness20.png", height = 800, width = 800)
ggplot(simulated_relatedness20, aes(x=sim)) +
  geom_histogram(color = "#ebefff", fill="#839ffe", alpha = 0.6, bins = 50) +
  geom_segment(aes(x = mean(spouses$kinship_20), 
                   xend = mean(spouses$kinship_20), 
                   y = 0, 
                   yend = Inf), 
               color = "#b77fff", 
               linetype = 'dashed', 
               linewidth = 1) + 
  #ggtitle("Simulated kinship whole sample") +
  labs(x = "Relatedness (20cM)", y = 'Frequency') +
  theme_classic() +
  theme(text = element_text(size = 30))
dev.off()



## PROBABILITY OF FINDING SPOUSAL VALUES 


# FOR IBD

# Set up parameters for simulation
set.seed(123)  # For reproducibility
n <- nrow(spousal_relatedness)
num_permutations <- 1000
permuted_means <- numeric(num_permutations)
observed_mean = mean(spousal_relatedness$sum[spousal_relatedness$type == 'spouse'])


# We randomly shuffle the type/relationship labels and calculate the mean relatedness of the "couples" in these permuted datasets. 
for (i in 1:num_permutations) {
  
  # Shuffle the relationship labels
  shuffled_relationship <- sample(spousal_relatedness$sum)
  
  # Calculate mean relatedness for the shuffled "couples"
  permuted_means[i] <- mean(shuffled_relationship[spousal_relatedness$type == 'spouse'])
}
permuted_means = data.frame(sim = permuted_means)

# Calculate mean and sd
mean_distribution = mean(permuted_means$sim) #2185.872
sd_distribution = sd(permuted_means$sim) #75.07019

# calculate Z-score
z_score <- (observed_mean - mean_distribution) / sd_distribution #-1.433838

# Calculate p-value for Z-score
p_value <- pnorm(z_score) #0.07580918

# Plot
png("/Users/inezd/PHD/Chapter_3/Plots/permuted_couples_IBD.png", width = 800, height = 800)
ggplot(permuted_means, aes(x=sim)) +
  geom_histogram(color = "#ebefff", fill="#839ffe", alpha = 0.6, bins = 50) +
  geom_segment(aes(x = mean(observed_mean), 
                   xend = mean(observed_mean), 
                   y = 0, 
                   yend = Inf), 
               color = "#b77fff", 
               linetype = 'dashed', 
               linewidth = 1) + 
  #ggtitle("Distribution of Permuted Mean Relatedness") +
  labs(x = "Mean sum IBD", y = 'Frequency') +
  theme_classic() +
  theme(text = element_text(size = 30))
dev.off()


# FOR KINSHIP 7

set.seed(123)  # For reproducibility
n <- nrow(spousal_relatedness)
num_permutations <- 1000
permuted_means <- numeric(num_permutations)
observed_mean = mean(spousal_relatedness$kinship_7[spousal_relatedness$type == 'spouse'])


# We randomly shuffle the type/relationship labels and calculate the mean relatedness of the "couples" in these permuted datasets. 
for (i in 1:num_permutations) {
  # Shuffle the relationship labels
  shuffled_relationship <- sample(spousal_relatedness$kinship_7)
  
  # Calculate mean relatedness for the shuffled "couples"
  permuted_means[i] <- mean(shuffled_relationship[spousal_relatedness$type == 'spouse'])
}
permuted_means = data.frame(sim = permuted_means)

# Calculate mean and sd
mean_distribution = mean(permuted_means$sim) #0.1379914
sd_distribution = sd(permuted_means$sim) #0.008358163

# calculate Z-score
z_score <- (observed_mean - mean_distribution) / sd_distribution #-1.6412

# Calculate p-value for Z-score
p_value <- pnorm(z_score) #0.05037798

# Save in separate dataframe
raute_spousalrel_per <- data.frame(dataset = 'spousal_rel', kinship = permuted_means$sim)

png("/Users/inezd/PHD/Chapter_3/Plots/permuted_couples_7cM.png", width = 800, height = 800)
ggplot(permuted_means, aes(x=sim)) +
  geom_histogram(color = "#ebefff", fill="#839ffe", alpha = 0.6, bins = 50) +
  geom_segment(aes(x = mean(observed_mean), 
                   xend = mean(observed_mean), 
                   y = 0, 
                   yend = Inf), 
               color = "#b77fff", 
               linetype = 'dashed', 
               linewidth = 1) + 
  #ggtitle("Distribution of Permuted Mean Relatedness") +
  labs(x = "Mean kinship", y = 'Frequency') +
  theme_classic() +
  theme(text = element_text(size = 30))
dev.off()




## WHAT ABOUT FOR THE WHOLE POPULATION?
population_relatedness <- merge(raute_IBD_summary, raute_7_kin[,c(3,7)], by = 'pair', all.x=T)
population_relatedness <- merge(population_relatedness, raute_20_kin[,c(3,7)], by = 'pair', all.x=T)
population_relatedness <- merge(population_relatedness, spouses[,c(3,13)], by = 'pair', all.x=T)
population_relatedness$type[is.na(population_relatedness$type)] <- 'non-spouse'

# Actually check: spouses versus non spouses
anova_pop <- aov(sum ~ type, data = population_relatedness)
summary(anova_pop) # 2.395  0.122
anova_pop <- aov(kinship_7 ~ type, data = population_relatedness)
summary(anova_pop) #3.028 0.0819


# Check: relatedness within and between patriclans
population_relatedness$id1 <- sub("\\-.*", "", population_relatedness$pair)
population_relatedness$id2 <- sub('.*-', '', population_relatedness$pair)
names(id_gen)[8] <- 'old_FID'
id_gen <- merge(id_gen, raute_fam[,c(2,8)], by = 'old_FID', all.x=T)
names(population_relatedness)[8] <- 'new_IID'
population_relatedness <- merge(population_relatedness, id_gen[,c(6,7,9)], by = 'new_IID', all.x=T)
names(population_relatedness)[c(1,9:11)] <- c('id1', 'new_IID', 'gottra_born.1', 'gottra_married.1')
population_relatedness <- merge(population_relatedness, id_gen[,c(6,7,9)], by = 'new_IID', all.x=T)
names(population_relatedness)[c(1,12:13)] <- c('id2', 'gottra_born.2', 'gottra_married.2')
population_relatedness$gottra_marry.1[is.na(population_relatedness$gottra_married.1)] <- population_relatedness$gottra_born.1[is.na(population_relatedness$gottra_married.1)]
population_relatedness$gottra_marry.1[!is.na(population_relatedness$gottra_married.1)] <- population_relatedness$gottra_married.1[!is.na(population_relatedness$gottra_married.1)]
population_relatedness$gottra_marry.2[is.na(population_relatedness$gottra_married.2)] <- population_relatedness$gottra_born.2[is.na(population_relatedness$gottra_married.2)]
population_relatedness$gottra_marry.2[!is.na(population_relatedness$gottra_married.2)] <- population_relatedness$gottra_married.2[!is.na(population_relatedness$gottra_married.2)]
population_relatedness$gottra_combo <- paste(pmin(population_relatedness$gottra_marry.1, population_relatedness$gottra_marry.2), pmax(population_relatedness$gottra_marry.1, population_relatedness$gottra_marry.2), sep = '-')
population_relatedness$gottra_check <- ifelse(population_relatedness$gottra_marry.1 == population_relatedness$gottra_marry.2, 'same', 'different')
gottra_relatedness <- population_relatedness %>% 
  group_by(gottra_check) %>% 
  summarise(mean_sum = mean(sum),
            mean_kin7 = mean(kinship_7),
            mean_kin20 = mean(kinship_20))
anova_got <- aov(sum ~ gottra_check, data = population_relatedness)
summary(anova_got) #16.06 6.2e-05 ***
glm_got <- glm(sum ~ gottra_check, data = population_relatedness)
summary(glm_got) # 4.008  6.2e-05 ***
anova_got <- aov(kinship_7 ~ gottra_check, data = population_relatedness)
summary(anova_got) #29.53 5.72e-08 ***

# perform bayesian model 
blm <- brm(kinship_7 ~ gottra_check,
            data = population_relatedness)
summary(blm)
plot(blm)
cond_effects_blm <- conditional_effects(blm, effects = "gottra_check")
cond_data_blm <- cond_effects_blm[[1]]  #
png("/Users/inezd/PHD/Chapter_3/Plots/brms_gottra.png", width = 800, height = 800)
ggplot(cond_data_blm, aes(x = gottra_check, y = estimate__, ymin = lower__, ymax = upper__, group = gottra_check)) +
  geom_pointrange(linewidth = 1) +  # Plot the point estimate with confidence intervals
  geom_point(aes(color = gottra_check),size = 5, alpha = 0.6) +  # Customize point appearance
  theme_classic() +
  scale_color_manual(values=c("#ffa635", "#839ffe")) +
  theme(text = element_text(size = 30),
        legend.position = "none" ) +
  labs(y = "Kinship coefficient", x = "Clan membership", color = 'Clan membership') +
  scale_y_continuous(limits = c(0.138,0.146), breaks = seq(0.138,0.146, by = 0.002))
dev.off()

# Are there distinct differences between the gottras?
blm2 <- brm(kinship_7 ~ gottra_combo,
           data = population_relatedness)
summary(blm2)
plot(blm2)
cond_effects_blm2 <- conditional_effects(blm2, effects = "gottra_combo")
cond_data_blm2 <- cond_effects_blm2[[1]]  # WOW! Subanshi is smaller and much more inbred
png("/Users/inezd/PHD/Chapter_3/Plots/brms_gottra2.png", width = 800, height = 800)
ggplot(cond_data_blm2, aes(x = gottra_combo, y = estimate__, ymin = lower__, ymax = upper__, group = gottra_combo)) +
  geom_pointrange(linewidth = 1) +  # Plot the point estimate with confidence intervals
  geom_point(aes(color = gottra_combo),size = 5, alpha = 0.6) +  # Customize point appearance
  theme_classic() +
  #scale_color_manual(values=c("#ffa635", "#839ffe")) +
  theme(text = element_text(size = 30),
        legend.position = "none" ) +
  labs(y = "Kinship coefficient", x = "Clan membership", color = 'Clan membership') +
  scale_y_continuous(limits = c(0.130,0.21), breaks = seq(0.130,0.21, by = 0.01))
dev.off()

# Are the differences still there if we ignore S-S 
blm3 <- brm(kinship_7 ~ gottra_check,
           data = population_relatedness[!(population_relatedness$gottra_combo == 'S-S'),])
summary(blm3)
plot(blm3)
cond_effects_blm3 <- conditional_effects(blm3, effects = "gottra_check")
cond_data_blm3 <- cond_effects_blm3[[1]]  #
png("/Users/inezd/PHD/Chapter_3/Plots/brms_gottra_checked.png", width = 800, height = 800)
ggplot(cond_data_blm3, aes(x = gottra_check, y = estimate__, ymin = lower__, ymax = upper__, group = gottra_check)) +
  geom_pointrange(linewidth = 1) +  # Plot the point estimate with confidence intervals
  geom_point(aes(color = gottra_check),size = 5, alpha = 0.6) +  # Customize point appearance
  theme_classic() +
  scale_color_manual(values=c("#ffa635", "#839ffe")) +
  theme(text = element_text(size = 30),
        legend.position = "none" ) +
  labs(y = "Kinship coefficient", x = "Clan membership", color = 'Clan membership') +
  scale_y_continuous(limits = c(0.138,0.146), breaks = seq(0.138,0.146, by = 0.002))
dev.off()

# What if we remove first degree relatedness?
blm4 <- brm(kinship_7 ~ gottra_check,
            data = population_relatedness[!(population_relatedness$kinship_7  >=0.177),])
summary(blm4)
plot(blm4)
cond_effects_blm4 <- conditional_effects(blm4, effects = "gottra_check")
cond_data_blm4 <- cond_effects_blm4[[1]]  #
png("/Users/inezd/PHD/Chapter_3/Plots/brms_gottra_nofirst.png", width = 800, height = 800)
ggplot(cond_data_blm4, aes(x = gottra_check, y = estimate__, ymin = lower__, ymax = upper__, group = gottra_check)) +
  geom_pointrange(linewidth = 1) +  # Plot the point estimate with confidence intervals
  geom_point(aes(color = gottra_check),size = 5, alpha = 0.6) +  # Customize point appearance
  theme_classic() +
  scale_color_manual(values=c("#ffa635", "#839ffe")) +
  theme(text = element_text(size = 30),
        legend.position = "none" ) +
  labs(y = "Kinship coefficient", x = "Clan membership", color = 'Clan membership') +
  scale_y_continuous(limits = c(0.131,0.136), breaks = seq(0.131,0.136, by = 0.001))
dev.off()


anova_got <- aov(kinship_20 ~ gottra_check, data = population_relatedness)
summary(anova_got) #35.74 2.38e-09 ***

observed_mean <- 0.124274  # Replace with your observed mean
t.test(population_relatedness$kinship_7, mu = observed_mean) #0.3078



# Start simulation
set.seed(123)  # For reproducibility
n <- nrow(population_relatedness)
num_permutations <- 1000
permuted_means <- numeric(num_permutations)
observed_mean = mean(population_relatedness$kinship_7[population_relatedness$type == 'spouse'])

# We randomly shuffle the type/relationship labels and calculate the mean relatedness of the "couples" in these permuted datasets. 
for (i in 1:num_permutations) {
  # Shuffle the relationship labels
  shuffled_relationship <- sample(population_relatedness$kinship_7)
  
  # Calculate mean relatedness for the shuffled "couples"
  permuted_means[i] <- mean(shuffled_relationship[population_relatedness$type == 'spouse'])
}
permuted_means = data.frame(sim = permuted_means)

# Calculate mean and sd
mean_distribution = mean(permuted_means$sim) # 0.1409953
sd_distribution = sd(permuted_means$sim) # 0.009649882

# calculate Z-score
z_score <- (observed_mean - mean_distribution) / sd_distribution # -1.732795

# Calculate p-value for Z-score
p_value <- pnorm(z_score) #0.04156606 

# Save in separate dataframe
raute_allrel_per <- data.frame(dataset = 'total_population', kinship = permuted_means$sim)

# Plot
png("/Users/inezd/PHD/Chapter_3/Plots/permuted_population_7cM.png", width = 800, height = 800)
ggplot(permuted_means, aes(x=sim)) +
  geom_histogram(color = "#fcd7e9", fill="#ff9dcc", alpha = 0.6, bins = 50) +
  geom_segment(aes(x = mean(observed_mean), 
                   xend = mean(observed_mean), 
                   y = 0, 
                   yend = Inf), 
               color = "#b77fff", 
               linetype = 'dashed', 
               linewidth = 1) + 
  #ggtitle("Distribution of Permuted Mean Relatedness") +
  labs(x = "Mean kinship", y = 'Frequency') +
  theme_classic() +
  theme(text = element_text(size = 30))
dev.off()



## try for IBD

set.seed(123)  # For reproducibility
n <- nrow(population_relatedness)
num_permutations <- 1000
permuted_means <- numeric(num_permutations)
observed_mean = mean(population_relatedness$sum[population_relatedness$type == 'spouse'])


# We randomly shuffle the type/relationship labels and calculate the mean relatedness of the "couples" in these permuted datasets. 
for (i in 1:num_permutations) {
  # Shuffle the relationship labels
  shuffled_relationship <- sample(population_relatedness$sum)
  
  # Calculate mean relatedness for the shuffled "couples"
  permuted_means[i] <- mean(shuffled_relationship[population_relatedness$type == 'spouse'])
}
permuted_means = data.frame(sim = permuted_means)

# Calculate mean and sd
mean_distribution = mean(permuted_means$sim) #2209.806
sd_distribution = sd(permuted_means$sim) #83.55252

# calculate Z-score
z_score <- (observed_mean - mean_distribution) / sd_distribution #-1.57473

# Calculate p-value for Z-score
p_value <- pnorm(z_score) #0.05765941

# Plot
png("/Users/inezd/PHD/Chapter_3/Plots/permuted_population_IBD.png", width = 800, height = 800)
ggplot(permuted_means, aes(x=sim)) +
  geom_histogram(color = "#fcd7e9", fill="#ff9dcc", alpha = 0.6, bins = 50) +
  geom_segment(aes(x = mean(observed_mean), 
                   xend = mean(observed_mean), 
                   y = 0, 
                   yend = Inf), 
               color = "#b77fff", 
               linetype = 'dashed', 
               linewidth = 1) + 
  #ggtitle("Distribution of Permuted Mean Relatedness") +
  labs(x = "Mean sum IBD", y = 'Frequency') +
  theme_classic() +
  theme(text = element_text(size = 30))
dev.off()



### NEXT STEP: ONLY THE POTENTIAL PARTNERS

# Not 1st degree relatedness, different-sex, above 16
raute_list <- raute_census_combined[(raute_census_combined$alive2022 == 'Yes'),c(1,7,11,13,14)]
write.csv(raute_list, "/Users/inezd/PHD/Genetic_Data/IBD_analyses/Analyses/raute_list.csv")
raute_pairs <- data.frame(t(combn(raute_list$id_2014,2)))
names(raute_list)[1] <- 'X1'
raute_pairs <- merge(raute_list, raute_pairs, by = 'X1')
names(raute_list)[1] <- 'X2'
raute_pairs <- merge(raute_list, raute_pairs, by = 'X2')
raute_pairs <- raute_pairs[!(raute_pairs$sex.x == raute_pairs$sex.y),]
raute_pairs <- raute_pairs[(raute_pairs$age2022.x > 16 & raute_pairs$age2022.y > 16),]

# Merge this with relatedness data
names(id_gen)[2] <- 'X1'
raute_pairs <- merge(raute_pairs, id_gen[,c(2,1)], by = 'X1')
names(id_gen)[2] <- 'X2'
raute_pairs <- merge(raute_pairs, id_gen[,c(2,1)], by = 'X2')
names(id_gen)[2] <- 'id_2014'
names(raute_fam)[1] <- 'ind.x'
raute_pairs <- merge(raute_pairs, raute_fam[,c(1,8)], by = 'ind.x', all.x=T)
names(raute_fam)[1] <- 'ind.y'
raute_pairs <- merge(raute_pairs, raute_fam[,c(1,8)], by = 'ind.y', all.x=T)
names(raute_fam)[1] <- 'V1'
raute_pairs <- raute_pairs[,c(13,14,5:12)]
raute_pairs$pair <- paste(pmin(raute_pairs$new_IID.x, raute_pairs$new_IID.y), pmax(raute_pairs$new_IID.x, raute_pairs$new_IID.y), sep = '-')
raute_pairs <- merge(raute_pairs, population_relatedness, by = 'pair', all.x=T)

# Add homozygosity
names(raute_pairs)[2] <- 'IID'
raute_pairs <- merge(raute_pairs, homozygosity[,c(2,6,7)], by = 'IID', all.x=T)
names(raute_pairs)[c(1,3)] <- c('new_IID.x', 'IID')
raute_pairs <- merge(raute_pairs, homozygosity[,c(2,6,7)], by = 'IID', all.x=T)
names(raute_pairs)[1] <- 'new_IID.y'
write.csv(raute_pairs, "/Users/inezd/PHD/Genetic_Data/IBD_analyses/Analyses/raute_pairs.csv")


## Allow only individuals with relatedness below 0.1768
raute_relatedness <- raute_pairs[(raute_pairs$kinship_7 < 0.1768),]

# Actually check: spouses versus non spouses
anova_eligible <- aov(sum ~ type, data = raute_relatedness)
summary(anova_eligible) # 1.508   0.22
anova_eligible <- aov(kinship_7 ~ type, data = raute_relatedness)
summary(anova_eligible) #3.025 0.0823

observed_mean <- 0.124274  # Replace with your observed mean
t.test(raute_relatedness$kinship_7, mu = observed_mean) #0.3078


#start random shuffling
set.seed(123)  # For reproducibility
n <- nrow(raute_relatedness)
num_permutations <- 1000
permuted_means <- numeric(num_permutations)
observed_mean = mean(raute_relatedness$sum[raute_relatedness$type == 'spouse'])


# We randomly shuffle the type/relationship labels and calculate the mean relatedness of the "couples" in these permuted datasets. 
for (i in 1:num_permutations) {
  # Shuffle the relationship labels
  shuffled_relationship <- sample(raute_relatedness$sum)
  
  # Calculate mean relatedness for the shuffled "couples"
  permuted_means[i] <- mean(shuffled_relationship[raute_relatedness$type == 'spouse'])
}
permuted_means = data.frame(sim = permuted_means)

# Calculate mean and sd
mean_distribution = mean(permuted_means$sim)
print(mean_distribution) #2138.475
sd_distribution = sd(permuted_means$sim)
print(sd_distribution) #51.97519

# calculate Z-score
z_score <- (observed_mean - mean_distribution) / sd_distribution
print(z_score) #-1.159058

# Calculate p-value for Z-score
p_value <- pnorm(z_score) 
print(p_value) #0.1232162

# For IBD
png("/Users/inezd/PHD/Chapter_3/Plots/permuted_possibility_IBD.png", width = 800, height = 800)
ggplot(permuted_means, aes(x=sim)) +
  geom_histogram(color = "#fcdeb8", fill="#ffa635", alpha = 0.6, bins = 50) +
  geom_segment(aes(x = mean(observed_mean), 
                   xend = mean(observed_mean), 
                   y = 0, 
                   yend = Inf), 
               color = "#b77fff", 
               linetype = 'dashed', 
               linewidth = 1) + 
  #ggtitle("Distribution of Permuted Mean Relatedness") +
  labs(x = "Mean sum IBD", y = 'Frequency') +
  theme_classic() +
  theme(text = element_text(size = 30))
dev.off()



# For kinship7
set.seed(123)  # For reproducibility
n <- nrow(raute_relatedness)
num_permutations <- 1000
permuted_means <- numeric(num_permutations)
observed_mean = mean(raute_relatedness$kinship_7[raute_relatedness$type == 'spouse'])


# We randomly shuffle the type/relationship labels and calculate the mean relatedness of the "couples" in these permuted datasets. 
for (i in 1:num_permutations) {
  # Shuffle the relationship labels
  shuffled_relationship <- sample(raute_relatedness$kinship_7)
  
  # Calculate mean relatedness for the shuffled "couples"
  permuted_means[i] <- mean(shuffled_relationship[raute_relatedness$type == 'spouse'])
}
permuted_means = data.frame(sim = permuted_means)

# Calculate mean and sd
mean_distribution = mean(permuted_means$sim)
print(mean_distribution) #0.1325096
sd_distribution = sd(permuted_means$sim)
print(sd_distribution) #0.004847791

# calculate Z-score
z_score <- (observed_mean - mean_distribution) / sd_distribution
print(z_score) #-1.698831

# Calculate p-value for Z-score
p_value <- pnorm(z_score) 
print(p_value) #0.04467551

# Save in separate dataframe
raute_elirel_per <- data.frame(dataset = 'eligible_pairs', kinship = permuted_means$sim)

# Plot
png("/Users/inezd/PHD/Chapter_3/Plots/permuted_possibility_7cm.png", width = 800, height = 800)
ggplot(permuted_means, aes(x=sim)) +
geom_histogram(color = "#fcdeb8", fill="#ffa635", alpha = 0.6, bins = 50) +
  geom_segment(aes(x = mean(observed_mean), 
                   xend = mean(observed_mean), 
                   y = 0, 
                   yend = Inf), 
               color = "#b77fff", 
               linetype = 'dashed', 
               linewidth = 1) + 
  #ggtitle("Distribution of Permuted Mean Relatedness") +
  labs(x = "Mean kinship", y = 'Frequency') +
  theme_classic() +
  theme(text = element_text(size = 30))
dev.off()


#### POPULATION COMPARISON ######

# First, for the means per three datasets
raute_spousalrel <- spousal_relatedness[,c(1,7,9)]
names(raute_spousalrel)[2] <- 'kinship'
raute_spousalrel$dataset <- 'spousal_rel'
raute_eligiblerel <- raute_relatedness[,c(1,15,17)]
names(raute_eligiblerel)[2] <- 'kinship'
raute_eligiblerel$dataset <- 'eligible_pairs'
raute_allrel <- population_relatedness[,c(3,7,9)]
names(raute_allrel)[2] <- 'kinship'
raute_allrel$dataset <- 'total_population'
raute_datasets <- rbind(raute_spousalrel, raute_eligiblerel, raute_allrel)
raute_datasets$population <- 'Raute'

# Check for a predictor?
spouses_kin <- raute_datasets[(raute_datasets$type == 'spouse'),c(1,2,4)]
spouses_kin$dataset <- 'spouses'
spouses_kin <- unique(spouses_kin)
spouses_kin <- rbind(spouses_kin, raute_datasets[(raute_datasets$type == 'non-spouse'),c(1,2,4)])
fit3 <- brm(kinship ~ dataset,
            data = spouses_kin,
            iter = 15000,          # Increase the number of iterations
            warmup = 5000,         # Increase the warmup period
            control = list(adapt_delta = 0.99))
fitla <- brm(kinship ~ dataset,
            data = spouses_kin,
            iter = 15000,          # Increase the number of iterations
            warmup = 5000,         # Increase the warmup period
            control = list(adapt_delta = 0.99),
            family = inverse.gaussian(link = "log"))
summary(fitla)
plot(conditional_effects(fitla, effects = 
                           "dataset"))

summary(fit3)
plot(conditional_effects(fit3, effects = "dataset"))
cond_effects2 <- conditional_effects(fit3, effects = "dataset")
cond_data2 <- cond_effects2[[1]]  #
ggplot(cond_data2, aes(x = dataset, y = estimate__, ymin = lower__, ymax = upper__, group = dataset)) +
  geom_pointrange(aes(color = dataset), linewidth = 1) +  # Plot the point estimate with confidence intervals
  geom_point(size = 5, alpha = 0.6) +  # Customize point appearance
  theme_classic() +
  scale_color_manual(values=c("#ffa635", "#839ffe", "#ff839f", "#cb8cff")) +
  theme(text = element_text(size = 25),) +
  labs(y = "Estimated FROH", x = "Population")

# Combine with Bayaka and Agta
rel_datasets <- rbind(agta_datasets, bayaka_datasets, raute_datasets)
rel_datasets$Relationship <- ifelse(rel_datasets$type == 'spouse', 'Couple', "Non-couple")
rel_datasets$dataset <- ifelse(rel_datasets$dataset == 'eligible_pairs', 'Eligible',
                               ifelse(rel_datasets$dataset == 'spousal_rel', 'Spousal', 'Total'))

# Expand observed_means to include all datasets for each population
observed_means <- expand.grid(
  dataset = c("Eligible", "Spousal", "Total"),
  population = factor(c("Agta", "BaYaka", "Raute"))
)

# Assign observed mean values for each population
observed_means$observed_mean <- rep(c(0.02088135, 0.005302667, 0.124274), each = 3)

# Plot
png("/Users/inezd/PHD/Chapter_3/Plots/BaYaka//all_pops_compare_plot.png", width =1600, height = 800)
ggplot(rel_datasets[(rel_datasets$Relationship == 'Non-couple'),], aes(x = dataset, y = kinship, group = dataset)) +
  geom_boxplot(position = position_dodge(width = 0.9), alpha = 0.6, aes(fill = dataset)) + 
  scale_fill_manual(values = c("Eligible" = "#ffa635", 
                               "Spousal" = "#839ffe",
                               "Total" = "#ff839f")) +  # Customize dataset colors
  theme_classic() +
  geom_segment(data = observed_means, aes(x = as.numeric(dataset) - 0.45, 
                                          xend = as.numeric(dataset) + 0.45, 
                                          y = observed_mean, 
                                          yend = observed_mean,
                                          group = population),  # Add 'population' to group aesthetics
               inherit.aes = FALSE,  # Prevent inheriting 'dataset' aesthetic from main plot
               linetype = "dashed", color = "purple", size = 1) +
  theme(text = element_text(size = 30),
        legend.position = 'none') +
  labs(x = "Dataset", y = "Kinship", fill = "Dataset") +
  facet_wrap(~ population, scales = "free_y") +
  scale_y_continuous(limits = c(0,0.35), breaks = seq(0,0.35, by = 0.05)) 
dev.off()


# Density plot
png("/Users/inezd/PHD/Chapter_3/Plots/BaYaka/all_pops_density_plot.png", width =1600, height = 800)
ggplot(rel_datasets, aes(x = kinship, fill = dataset)) +
  geom_density(alpha = 0.6) +  # Density plot for distribution
  facet_wrap(~ population) +  # Facet by population
  geom_vline(data = observed_means, aes(xintercept = observed_mean), linetype = "dashed", color = "#b77fff", size = 1) +
  scale_fill_manual(values = c("#ffa635", "#839ffe", "#ff839f")) +  # Customize colors for datasets
  theme_classic() +
  labs(x = "Kinship", y = "Density", fill = "Dataset") +
  theme(text = element_text(size = 30))
dev.off()

## Then for permutations

# Combine raute files
raute_permutations <- rbind(raute_elirel_per, raute_allrel_per, raute_spousalrel_per)
raute_permutations$population <- 'Raute'
agta_permutations$population <- 'Agta'
bayaka_permutations$population <- 'BaYaka'
all_permutations <- rbind (raute_permutations, agta_permutations, bayaka_permutations)
all_permutations$dataset <- ifelse(all_permutations$dataset == 'eligible_pairs', 'Eligible',
                               ifelse(all_permutations$dataset == 'spousal_rel', 'Spousal', 'Total'))

# Model
fit5 <- brm(kinship ~ dataset + (1|population),
            data = all_permutations,
            iter = 10000)
summary(fit4)
plot(conditional_effects(fit4, effects = "dataset:population"))
cond_effects4 <- conditional_effects(fit3, effects = "dataset:population")
cond_data4 <- cond_effects4[[1]]  #

# Create a data frame for the horizontal lines
line_data <- data.frame(
  population = c("Agta", "BaYaka", "Raute"),
  y_value = c(0.0209, 0.00530, 0.124)
)
png("/Users/inezd/PHD/Chapter_3/Plots/BaYaka/model_permutations_plot.png", width =1600, height = 800)
ggplot(cond_data4, aes(x = dataset, y = estimate__, ymin = lower__, ymax = upper__, group = dataset)) +
  geom_pointrange(aes(color = dataset), linewidth = 1) +  # Plot the point estimate with confidence intervals
  geom_point(size = 5, alpha = 0.6, aes(color = dataset)) +  # Customize point appearance
  theme_classic() +
  scale_color_manual(values=c("#ffa635", "#839ffe", "#ff839f")) +
  theme(text = element_text(size = 30),
        legend.position = 'none') +
  labs(y = "Kinship", x = "Dataset", color = 'Dataset') +
  geom_hline(data = line_data, aes(yintercept = y_value), 
             linetype = "dashed", color = "purple", size = 1) +  # Add horizontal lines per population
  facet_wrap(~ population, scales = "free_y")  # Facet by population with free y-scales

dev.off()


# Density plot
png("/Users/inezd/PHD/Chapter_3/Plots/BaYaka//density_permutations_plot.png", width =1600, height = 800)
ggplot(all_permutations, aes(x = kinship, fill = dataset)) +
  geom_density(alpha = 0.6) +  # Density plot for distribution
  facet_wrap(~ population) +  # Facet by population
  geom_vline(data = observed_means, aes(xintercept = observed_mean), linetype = "dashed", color = "#b77fff", size = 1) +
  scale_fill_manual(values = c("#ffa635", "#839ffe", "#ff839f")) +  # Customize colors for datasets
  theme_classic() +
  labs(x = "Kinship", y = "Density", fill = "Dataset") +
  theme(text = element_text(size = 30),
        legend.position = 'none')
dev.off()
