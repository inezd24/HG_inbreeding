# Script developed by Inez Derkx, summer 2024
# Script aims to examine list of spouses and their relatedness alongside relatedness between all those 22 individuals. 

##################################### IMPORT LIBRARIES #####################################

library(ggplot2)
library(reshape2)
library(readr)
library(dplyr)

##################################### IMPORT DATA #####################################

spouses <- read.csv("/Users/inezd/PHD/Chapter_3/R_Files/spousal_information.csv", header = T, row.names = 1) # Df with info of all spouses
raute_IBD_summary <- read.csv("/Users/inezd/PHD/Chapter_3/R_Files/raute_IBD_summary.csv", header = T, row.names = 1) # Df with IBD-based relatedness
raute_7_kin <- read.delim("~/PHD/Genetic_Data/IBIS/raute_maf_7.coef",
                          col.names = c("sample_1", "sample_2", "kinship_7", "IBD2_7", "segment_count", "degree_of_relatedness")) # Df with kinship coefficients
id_gen # Df with different IDs of individuals
raute_fam # Df with different IDs of individuals
homozygosity # Df with FROH measures

##################################### TRANSFORM DATA #####################################

## AMONG INDIVIDUALS WHO ARE SPOUSES

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

# Examine the probability of finding the observed values by permuting relatedness (see function)
permuted_means <- permute_relatedness(spousal_relatedness)

# Save in separate dataframe for later merging
raute_spousalrel_per <- data.frame(dataset = 'spousal_rel', kinship = permuted_means$permuted_means_df$sim)

# Plot the results
png("/Users/inezd/PHD/Chapter_3/Plots/permuted_couples_7cM.png", width = 800, height = 800)
ggplot(permuted_means$permuted_means_df, aes(x=sim)) +
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

## AMONG TOTAL POPULATION

# What does relatedness look like in the entire population
population_relatedness <- raute_IBD_summary %>% 
  left_join(raute_7_kin %>% select(kinship = kinship_7, pair), by = 'pair') %>%
  left_join(spouses %>% select(pair, type), by = 'pair') %>%
  mutate(type = replace_na(type, "non-spouse")) %>%
  na.omit() %>%
  mutate(id1 = sub("\\-.*", "", pair),
         id2 = sub(".*-", "", pair)) %>%
  left_join(id_gen %>% select(id1 = new_IID, gottra_born.1 = gottra_born, gottra_married.1 = gottra_married), by = 'id1') %>%
  left_join(id_gen %>% select(id2 = new_IID, gottra_born.2 = gottra_born, gottra_married.2 = gottra_married), by = 'id2') %>%
  mutate(gottra_marry.1 = if_else(is.na(gottra_married.1), gottra_born.1, gottra_married.1),
         gottra_marry.2 = if_else(is.na(gottra_married.2), gottra_born.2, gottra_married.2),
         gottra_combo = paste(pmin(gottra_marry.1, gottra_marry.2), 
                              pmax(gottra_marry.1, gottra_marry.2), sep = "-"),
         gottra_check = if_else(gottra_marry.1 == gottra_marry.2, "same", "different"))

# Visualize observed true mean versus distribution between spouses
png("/Users/inezd/PHD/Chapter_3/Plots/oopulation_kin_avg.png", height = 800, width = 800)
population_kin <- ggplot(population_relatedness, aes(x = kinship)) +
  geom_histogram(color = "#ebefff", fill="#839ffe", alpha = 0.6, bins = 50) +
  geom_density(alpha=.2, fill="#FF6666") +
  geom_segment(aes(x = 0.12, xend = 0.12, y = 0, yend = Inf), color = "#b77fff", linetype = 'dashed', linewidth = 1) + # amongst spouses // purple
  theme_classic() +
  theme(text = element_text(size = 30)) +
  scale_x_continuous(limits = c(0.05, 0.30), breaks = seq(0.05, 0.30, by = 0.05)) +
  labs(title = "", x = 'Kinship coefficient', y = 'Frequency')
spousal_kin
dev.off()

# Examine the probability of finding the observed values by permuting relatedness (see function)
permuted_means <- permute_relatedness(population_relatedness)

# Save in separate dataframe for later merging
raute_allrel_per <- data.frame(dataset = 'total_population', kinship = permuted_means$permuted_means_df$sim)

# Plot the results
png("/Users/inezd/PHD/Chapter_3/Plots/permuted_population.png", width = 800, height = 800)
ggplot(permuted_means$permuted_means_df, aes(x=sim)) +
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

## AMONG ELIGIBLE INDIVIDUALS / POTENTIAL PARTNERS

# Create new dataframe and save for later use
raute_list <- raute_census_combined %>%
  filter(alive2022 == "Yes") %>%
  select(id_2014, sex, age2022, gottra_born, gottra_married) 
write.csv(raute_list, "/Users/inezd/PHD/Genetic_Data/IBD_analyses/Analyses/raute_list.csv")

# From this dataframe, generate pairs
raute_pairs <- combn(raute_list$id_2014, 2) %>%
  t() %>%
  as.data.frame() %>%
  rename(X1 = V1, X2 = V2) %>%
  left_join(raute_list %>% 
              select(X1 = id_2014, sex_x = sex, age_x = age2022), by = 'X1') %>%
  left_join(raute_list %>% 
              select(X2 = id_2014, sex_y = sex, age_y = age2022), by = 'X2') %>%
  filter(sex_x != sex_y, age_x > 16, age_y > 16) %>% # Only individuals of opposite sex and older than 16 
  left_join(id_gen %>% select(X1 = id_2014, ind.x= ind), by = 'X1') %>%
  left_join(id_gen %>% select(X2 = id_2014, ind.y= ind), by = 'X2') %>%
  left_join(raute_fam %>% select(ind.x = V1, id1 = new_IID), by = 'ind.x') %>%
  left_join(raute_fam %>% select(ind.y = V1, id2 = new_IID), by = 'ind.y') %>%
  select(c(9,10,3:6)) %>%
  mutate(pair = paste(pmin(id1, id2), pmax(id1, id2), sep = '-')) %>%
  left_join(population_relatedness, by = c('id1', 'id2', 'pair')) %>% # add relatedness values and gottra information
  left_join(homozygosity %>% select(id1 = IID, FROH.x = FROH), by = 'id1') %>%
  left_join(homozygosity %>% select(id2 = IID, FROH.y = FROH), by = 'id2') %>%
  filter(kinship < 0.1768) #496 dyads
write.csv(raute_pairs, "/Users/inezd/PHD/Genetic_Data/IBD_analyses/Analyses/raute_pairs.csv")

# Visualize observed true mean versus distribution between spouses
png("/Users/inezd/PHD/Chapter_3/Plots/eligible_kin_avg.png", height = 800, width = 800)
eligible_kin <- ggplot(raute_pairs, aes(x = kinship)) +
  geom_histogram(color = "#ebefff", fill="#839ffe", alpha = 0.6, bins = 50) +
  geom_density(alpha=.2, fill="#FF6666") +
  geom_segment(aes(x = 0.12, xend = 0.12, y = 0, yend = Inf), color = "#b77fff", linetype = 'dashed', linewidth = 1) + # amongst spouses // purple
  theme_classic() +
  theme(text = element_text(size = 30)) +
  scale_x_continuous(limits = c(0.05, 0.30), breaks = seq(0.05, 0.30, by = 0.05)) +
  labs(title = "", x = 'Kinship coefficient', y = 'Frequency')
spousal_kin
dev.off()

# Examine the probability of finding the observed values by permuting relatedness (see function)
permuted_means <- permute_relatedness(raute_pairs)

# Save in separate dataframe
raute_elirel_per <- data.frame(dataset = 'eligible_pairs', kinship = permuted_means$permuted_means_df$sim)

# Plot
png("/Users/inezd/PHD/Chapter_3/Plots/permuted_possibility.png", width = 800, height = 800)
ggplot(permuted_means$permuted_means_df, aes(x=sim)) +
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

## End of Script

