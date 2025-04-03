# Script developed by Inez Derkx
# Script aims to examine list of spouses and their relatedness alongside relatedness between all those 22 individuals. 

##################################### IMPORT LIBRARIES #####################################

library(ggplot2)
library(reshape2)
library(readr)
library(dplyr)

##################################### IMPORT DATA #####################################

# Import data
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


##################################### CREATE COLOUR LIST #####################################

# Add colour list to my colour palette
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


##################################### SPOUSES FILE #####################################

## STEP 1: GET RELATEDNESS FOR SPOUSES

# Create a file with the relevant ids for spouses
spouses <- raute_genealogies_complete %>%
  filter(!is.na(partner_id), partner_alive == 'yes', alive == 'yes') %>%
  select(id_2014_1 = id_2014, id = old_id, name = name, partner_name = partner_name, old_id = partner_id) %>%
  left_join(raute_genealogies_complete %>% select(old_id = old_id, id_2014_2 = id_2014), by = "old_id") %>%
  select(name = name, partner_name = partner_name, id_2014_1 = id_2014_1, id_2014_2 = id_2014_2) %>%
  arrange(name) %>%
  mutate(
    pair = paste(pmin(id_2014_1, id_2014_2), pmax(id_2014_1, id_2014_2), sep = '-'),
    id1 = pmin(id_2014_1, id_2014_2),
    id2 = pmax(id_2014_1, id_2014_2)
  ) %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(id1, id2) %>%
  arrange(id1, id2) %>%
  mutate_at(c('id1', 'id2'), as.numeric) %>%
  left_join(id_gen %>% select(ind1 = ind, id1 = id_2014), by = "id1") %>%
  left_join(id_gen %>% select(ind2 = ind, id2 = id_2014), by = "id2") %>%
  na.omit() %>%
  left_join(raute_fam %>% select(ind1 = V1, sample_1 = new_IID), by = 'ind1') %>% # if you get a warning: make sure to ungroup:
  left_join(raute_fam %>% select(ind2 = V1, sample_2 = new_IID), by = 'ind2') %>% # raute_fam <- raute_fam %>% ungroup(V3)
  select(3,4,1,2,5,6) %>%
  mutate(pair = paste(pmin(sample_1, sample_2), pmax(sample_1, sample_2), sep = '-'))

# Add IBD sharing
raute_IBD_summary <- raute_only_IBD %>%
  mutate(pair = paste(pmin(sample_1, sample_2), pmax(sample_1, sample_2), sep = '-')) %>%
  group_by(pair) %>%
  summarise(
    n = n(),
    sum = sum(length_cM),
    mean = mean(length_cM),
    .groups = "drop")

# Summary stats
raute_IBD_summary %>% summarise(
  mean_sum = round(mean(sum),3), # Mean sum: 2207.67
  mean_mean = round(mean(mean),3), # Mean mean: 14.91
  mean_n = round(mean(n),3)) # Mean number: 169.17

# Merge this with the spousal info
spouses <- spouses %>% 
  inner_join(raute_IBD_summary, by = "pair")

# Summary stats
spouses %>% summarise(
  mean_sum = round(mean(sum),3), # Mean sum: 2078.2
  mean_mean = round(mean(mean),3), # Mean mean: 11.5
  mean_n = round(mean(n),3)) # Mean number: 182.64

# Add kinship coefficient
raute_7_kin <- raute_7_kin %>%
  mutate(pair = paste(pmin(sample_1, sample_2), pmax(sample_1, sample_2), sep = '-'))
spouses <- spouses %>% 
  inner_join(raute_7_kin %>% select(pair = pair, kinship = kinship_7), by = "pair") %>%
  mutate(type = 'spouse')
mean(spouses$kinship) # Mean is 0.124274

# Create function for expected inbreeding coefficient of offspring
expected_inbreeding_offspring <- function(father, mother, kinship) {
  return(0.5 * (father + mother) + 0.5 * kinship)
}

# Finally, add to dataframe
spouses <- spouses %>%
  inner_join(roh_compare %>% select(sample_1 = IID, FROH_1 = FROH), by = "sample_1") %>%# Only select columns 'IID' and 'FROH'
  inner_join(roh_compare %>% select(sample_2 = IID, FROH_2 = FROH), by = "sample_2") %>%
  mutate(inbreeding_coef = 0.5*(FROH_1 + FROH_2) + (0.5*kinship))
mean(spouses$inbreeding_coef) # Mean is 0.292846

# Save file
write.csv(spouses, "/Users/inezd/Documents/Science/Raute/Chapter_3/R_Files/spousal_information.csv")
write.csv(raute_IBD_summary, "/Users/inezd/Documents/Science/Raute/R_Files/raute_IBD_summary.csv")

# Supplementary Figure 4
spouses$couple <- factor(spouses$pair, levels = spouses$pair[order(spouses$kinship)])
png("/Users/inezd/Documents/Science/Raute/Chapter_3/Plots/compare_spouses.png", height = 1000, width = 1000)
ggplot(spouses, aes(x = kinship, y = couple)) +
  geom_point(size = 5, color = "#ff839f", alpha = 0.7) +  # Customize point size and color
  geom_vline(xintercept = 0.141, color = '#b77fff', linetype="dotted", size = 2) +
  theme_classic() +
  scale_x_continuous(limits = c(0.10, 0.15)) +
  labs(x = "Relatedness", y = "Couples") +
  theme(axis.text.y = element_blank(),
        text = element_text(size = 40)) 
dev.off()
spouses$couple <- NULL


##################################### ELIGIBLE POPULATION #####################################

# Create datafile with eligible pairs of Raute individuals
raute_list <- raute_census_combined %>%
  filter(alive2022 == "Yes") %>%
  select(id_2014, sex, age2022, gottra_born, gottra_married) 
write.csv(raute_list, "/Users/inezd/Documents/Science/Raute/Chapter_3/R_Files/raute_list.csv")

# From this dataframe, generate pairs
raute_eligible <- combn(raute_list$id_2014, 2) %>%
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
write.csv(raute_eligible, "/Users/inezd/Documents/Science/Raute/Chapter_3/R_Files/raute_eligible.csv")


##################################### TOTAL POPULATION #####################################

# Total population relatedness
raute_total <- raute_IBD_summary %>% 
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
write.csv(raute_total, "/Users/inezd/Documents/Science/Raute/Chapter_3/R_Files/raute_total.csv")


##################################### PERMUTE RELATEDNESS #####################################

# Use function from script 'Function_Permute_relatedness.R'

# Apply function
permuted_means <- permute_relatedness(raute_total)

# Save in separate dataframe for later merging
raute_permuted <- data.frame(Population = 'Raute', kinship = permuted_means$permuted_means_df$sim)
write.csv(raute_permuted, "/Users/inezd/Documents/Science/Raute/Chapter_3/R_Files/raute_permuted.csv")



## End of Script




















