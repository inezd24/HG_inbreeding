# Script developed by Inez Derkx, summer 2024
# Script aims to examine list of spouses and their relatedness alongside relatedness between all those 22 individuals. 

##################################### IMPORT LIBRARIES #####################################

library(ggplot2)
library(reshape2)
library(readr)
library(dplyr)

##################################### IMPORT DATA #####################################

spouses <- read.csv("/Users/inezd/PHD/Chapter_3/R_Files/spousal_information.csv")

##################################### TRANSFORM DATA #####################################

# Create list of single individuals who are spouses
spouse_ind <- spouses %>%
  select(pair, sample_1, sample_2) %>%
  pivot_longer(cols = c(sample_1, sample_2), values_to = "sample_1", names_to = "variable") %>%
  distinct(sample_1)

# Create pairwise list of combinations of individuals who are spouses
spousal_relatedness <- expand_grid(sample_1 = spouse_ind$sample_1, sample_2 = spouse_ind$sample_1) %>%
  mutate(id1 = pmin(sample_1, sample_2),
         id2 = pmax(sample_1, sample_2),
         pair = paste(pmin(id1, id2), pmax(id1, id2), sep = '-')) %>%
  select(X1 = id1, X2 = id2, pair) %>%
  distinct() %>%
  left_join(raute_IBD_summary, by = 'pair') %>%
  left_join(raute_7_kin %>% select(kinship = kinship_7, pair), by = 'pair') %>%
  left_join(spouses %>% select(pair, type), by = 'pair') %>%
  mutate(type = replace_na(type, "non-spouse")) %>%
  na.omit() # 11 spouses, 220 non-spouses

# Check summary statistics (first alone, then for spouses vs. non-spouses)
spousal_relatedness %>% 
  summarise(
    mean_sum = mean(sum), # 2184.90
    mean_kin = mean(kinship)) # 0.14
spousal_relatedness %>%   
  group_by(type) %>% 
  summarise(mean_IBD = mean(sum),
                               n_IBD = mean(n),
                               mean_kin = mean(kinship)) 
observed_mean <- 0.124274 # Mean relatedness of true couples

# Visualize observed true mean versus distribution between spouses
png("/Users/inezd/PHD/Chapter_3/Plots/spousal_kin_avg.png", height = 800, width = 800)
spousal_kin <- ggplot(spousal_relatedness, aes(x = kinship)) +
  geom_histogram(color = "#ebefff", fill="#839ffe", alpha = 0.6, bins = 50) +
  geom_density(alpha=.2, fill="#FF6666") +
  geom_segment(aes(x = 0.12, xend = 0.12, y = 0, yend = Inf), color = "#b77fff", linetype = 'dashed', linewidth = 1) + # amongst spouses // purple
  theme_classic() +
  theme(text = element_text(size = 30)) +
  scale_x_continuous(limits = c(0.05, 0.30), breaks = seq(0.05, 0.30, by = 0.05)) +
  labs(title = "", x = 'Kinship coefficient', y = 'Frequency')
spousal_kin
dev.off()

# Let's go a bit deeper: simulate relatedness within this sample
set.seed(123) # For reproducibility
simulated_relatedness <- replicate(10000, { # Perform 10000 permutations
  sampled_pair <- spousal_relatedness[sample(nrow(spousal_relatedness), 11), ] # Sample 11 individuals (same number as there are spousal pairs)
  return(mean(sampled_pair[["kinship"]])) # Return the mean relatedness between them
})

# Obtain the mean and confidence intervals
simulated_relatedness <- data.frame(sim = simulated_relatedness)
mean_sim <- mean(simulated_relatedness$sim) #mean 
ci_sim <- quantile(simulated_relatedness$sim, c(0.025, 0.975)) #CI

# Plot distribution of simulated relatedness
png("/Users/inezd/PHD/Chapter_3/Plots/simulated_relatedness.png", height = 800, width = 800)
ggplot(simulated_relatedness, aes(x = sim)) +
  geom_histogram(color = "#ebefff", fill = "#839ffe", alpha = 0.6, bins = 50) +
  geom_vline(aes(xintercept = mean(spouses$kinship)), 
             color = "#b77fff", linetype = "dashed", linewidth = 1) + 
  geom_vline(aes(xintercept = ci_sim[1]), color = "red", linetype = "dotted") +
  geom_vline(aes(xintercept = ci_sim[2]), color = "red", linetype = "dotted") +
  labs(x = "Relatedness", y = "Frequency") +
  theme_classic() +
  theme(text = element_text(size = 20))
dev.off()


---------


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
