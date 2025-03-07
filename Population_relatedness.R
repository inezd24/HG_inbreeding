# Script developed by Inez Derkx, summer 2024
# Script aims to examine list of all individuals in a population for comparison. 

##################################### IMPORT LIBRARIES #####################################

library(ggplot2)
library(reshape2)
library(readr)
library(dplyr)

##################################### IMPORT DATA #####################################

raute_IBD_summary <- read.csv(raute_IBD_summary, "/Users/inezd/PHD/Chapter_3/R_Files/raute_IBD_summary.csv", header = T, sep = ',')
raute_7_kin <- read.delim("~/PHD/Genetic_Data/IBIS/raute_maf_7.coef",
                          col.names = c("sample_1", "sample_2", "kinship_7", "IBD2_7", "segment_count", "degree_of_relatedness"))

##################################### TRANSFORM DATA #####################################

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
















