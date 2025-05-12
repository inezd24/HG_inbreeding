# Script developed by Inez Derkx, summer 2024 and edited April 2025
# Script aims to examine population relatedness in the Agta, Bayaka, and Raute
# This script requires a function outlined in the script 'Function_Permute_relatedness.R' as well as data from 
# the three populations, which is transformed in the scripts 'Bayaka_Agta_datasets.R' and 'Raute_datasets.R'.

##################################### IMPORT LIBRARIES #####################################

library(ggplot2)
library(reshape2)
library(readr)
library(dplyr)

##################################### IMPORT DATA #####################################

bayaka_permuted <- read.csv("/Users/inezd/Documents/Science/Raute/Chapter_3/R_Files/bayaka_permuted.csv", header = T, row.names = 1)
agta_permuted <- read.csv("/Users/inezd/Documents/Science/Raute/Chapter_3/R_Files/agta_permuted.csv", header = T, row.names = 1)
raute_permuted <- read.csv("/Users/inezd/Documents/Science/Raute/Chapter_3/R_Files/raute_permuted.csv", header = T, row.names = 1)

agta_total <- read.csv("/Users/inezd/Documents/Science/Raute/Chapter_3/R_Files/agta_total.csv", header = T, row.names = 1)
bayaka_total <- read.csv("/Users/inezd/Documents/Science/Raute/Chapter_3/R_Files/bayaka_total.csv", header = T, row.names = 1)
raute_total <- read.csv("/Users/inezd/Documents/Science/Raute/Chapter_3/R_Files/raute_total.csv", header = T, row.names = 1)


##################################### FIGURE 2A #####################################

# Save the means
observed_means <- data.frame(
  population = c("Agta", "BaYaka", "Raute"),
  observed_mean = c(0.02088135, 0.005302667, 0.124274), # this is with the Agta outlier
  observed_mean2 <- c(0.01749632, 0.005302667, 0.124274)) # This is without the Agta outlier

# Total population distribution
total_datasets <- rbind(agta_total, bayaka_total, raute_total)
png("/Users/inezd/Documents/Science/Raute/Chapter_3/Plots/pops_compare_total.png", height = 1000, width = 1000)
total_rel_compare <- ggplot(rel_datasets_total, aes(x = kinship, fill = population)) +
  geom_density(color = "#ebefff", aes(fill=population), alpha = 0.6) +
  geom_vline(data = observed_means, aes(xintercept = observed_mean2, color = population), 
             linetype = "dashed", size = 3) +
  scale_fill_manual(values = c("Agta" = "#ffa635", 
                               "BaYaka" = "#839ffe",
                               "Raute" = "#ff839f")) +
  scale_color_manual(values = c("Agta" = "#ffa635", 
                               "BaYaka" = "#839ffe",
                               "Raute" = "#ff839f")) +
  theme_classic() +
  labs(title = "", x = 'Kinship coefficient', y = 'Frequency') +
  theme(text = element_text(size = 40),
        legend.position = 'none')
total_rel_compare
dev.off()


##################################### FIGURE 2B #####################################

# Combine files
eligible_permutations <- rbind(agta_permuted, bayaka_permuted, raute_permuted)

png("/Users/inezd/PHD/Chapter_3/Plots/pops_compare_eligible.png", height = 1000, width = 1000)
ggplot(eligible_permutations, aes(x = population, y = kinship, group = population)) +
  geom_boxplot(position = position_dodge(width = 0.9), alpha = 0.6, aes(fill = population)) + 
  scale_fill_manual(values = c("Agta" = "#ffa635", 
                               "BaYaka" = "#839ffe",
                               "Raute" = "#ff839f")) +  # Customize dataset colors
  geom_segment(data = observed_means, aes(x = as.numeric(population) - 0.45, 
                                          xend = as.numeric(population) + 0.45, 
                                          y = observed_mean2, 
                                          yend = observed_mean2,
                                          group = population,
                                          color = population),  # Add 'population' to group aesthetics
               inherit.aes = FALSE,  # Prevent inheriting 'dataset' aesthetic from main plot
               linetype = "dashed", size = 1) +
  scale_color_manual(values = c("Agta" = "#ffa635", 
                                "BaYaka" = "#839ffe",
                                "Raute" = "#ff839f")) +
  theme_classic() +
  theme(text = element_text(size = 40),
        legend.position = 'none') +
  labs(x = "Population", y = "Kinship", fill = "Population") +
  scale_y_continuous(limits = c(0,0.16), breaks = seq(0,0.16, by = 0.05)) 
dev.off()
