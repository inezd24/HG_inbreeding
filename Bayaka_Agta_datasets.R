# Script to edit the 'raw' BaYaka data (collected between 2013-2014), edited April 2025

##################################### IMPORT LIBRARIES #####################################

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

##################################### IMPORT DATA #####################################

bayaka_data <- read_excel("/Users/inezd/Documents/Science/Raute/Data/Bayaka_2018/Congo_2018_JULY2021.xlsx") # BaYaka general information file by individual
agta_fam <- read.table("/Users/inezd/Documents/Science/Raute/Genetic_Data/fam_files/AgtaHunterGatherer.fam") # Agta and BaYaka ID dataframe linking genetic and db ids
bayaka_coef <- read.delim("/Users/inezd/Documents/Science/Raute/Genetic_Data/IBIS/BaYaka.coef", 
                          col.names = c("sample_1", "sample_2", "kinship", "IBD2", "segment_count", "degree_of_relatedness")) # BaYaka output from IBIS (https://github.com/williamslab/ibis)
agta_data <- read_excel("/Users/inezd/Documents/Science/Raute/Genetic_Data/data_others/agta_data/Agta_Individuals2024-05-30_07_36_18.xlsx") # Agta general information file by individual
agta_new_ids <- read_excel("/Users/inezd/Documents/Science/Raute/Genetic_Data/data_others/agta_data/agta_person.xlsx") # extra info about ids from Agta individuals
agta_spouse <- read_excel("/Users/inezd/Documents/Science/Raute/Genetic_Data/data_others/agta_data/agta_spouses.xlsx") # Information about spouses
agta_coef <- read.delim("/Users/inezd/Documents/Science/Raute/Genetic_Data/IBIS/Agta.coef", 
                           col.names = c("sample_1", "sample_2", "kinship", "IBD2", "segment_count", "degree_of_relatedness")) # Agta output from IBIS

##################################### SPOUSES FILES #####################################

# Change ID dataframe
AGTA_BAYAKA <- agta_fam %>%
  group_by(V1) %>%
  mutate(new_IID = paste0(V1, row_number())) %>%
  ungroup() %>%
  mutate(new_pheno = 0) %>%
  rename(old_FID = V2, parent1 = V3, parent2 = V4, sex = V5, pheno = V6) %>%
  mutate(V3 = V1) %>%
  select(V1, old_FID, parent1, parent2, sex, pheno, V3, new_IID, new_pheno)
write.csv(AGTA_BAYAKA, "/Users/inezd/Documents/Science/Raute/Chapter_3/R_Files/AGTA_BAYAKA.csv")

# Create genealogies file for Agta individuals
agta_genealogies <- agta_data %>%
select(short_ID = ID, 'ID of Mother', 'ID of Father') %>% # select only necessary columns: person's id, mother id and father id. 
left_join(agta_new %>% select(person_id, short_ID = id_old), by = 'short_ID') %>% # join with the id file to add the new id
left_join(agta_spouse %>% select(person_id = pers_id, spouse_id_old), by = 'person_id') %>% # join with file with info on who is whose spouse
select(short_ID, old_ID = person_id, mother_id = 'ID of Mother', father_id = 'ID of Father', spouse_id_old) %>% # rename them
mutate(spouse_id_old = ifelse(grepl("F", spouse_id_old), NA, spouse_id_old)) %>% # delete non-relevant spouse ids
rename(spouse_id = spouse_id_old) %>%
distinct()
write.csv(agta_genealogies, "/Users/inezd/Documents/Science/Raute/Chapter_3/R_Files/agta_genealogies.csv")

# Add pairs to relatedness files (both Agta and BaYaka)
add_pair <- function(coef_file) {
  library(dplyr)
  coef_file2 <- coef_file %>%
    mutate(pair = paste(pmin(sample_1, sample_2), 
                        pmax(sample_1, sample_2), 
                        sep = '-'))
  return(coef_file2)
}
bayaka_coef <- add_pair(bayaka_coef)
agta_coef <- add_pair(agta_coef)

# Create spouses file
create_spouse_file <- function(data_file, id_file, coef_file){ # the data file needs to have an id column called 'short_ID' and a spouse column called 'spouse_id'
  out_file <- data_file %>%
    select(short_ID, spouse_id) %>% # select columns with ID of individual and ID of their spouse
    filter(!is.na(spouse_id)) %>% # filter out rows where there is no spouse
    left_join(id_file %>% select(short_ID = old_FID, new_IID), by = 'short_ID') %>% # add new id of individual
    left_join(id_file %>% select(spouse_id = old_FID, new_IID), by = 'spouse_id') %>% # add new id of spouse
    rename(person_id = short_ID, sample_1 = new_IID.x, sample_2 = new_IID.y) %>% # rename them
    mutate(id1 = pmin(sample_1, sample_2), # create alphabetic id
           id2 = pmax(sample_1, sample_2), # create alaphabetic id
           pair = paste(pmin(id1, id2), pmax(id1, id2), sep = '-')) %>% # create pair
    inner_join(coef_file %>% select(kinship, pair), by = 'pair') %>% # join with relatedness file by pair
    mutate(type = 'spouse') %>% # add variable that shows these are spouses
    select(pair, kinship, type) %>%
    distinct() 
  return(out_file) # return output file
}

# Apply function
bayaka_spouses <- create_spouse_file(bayaka_data, AGTA_BAYAKA, bayaka_coef)
agta_spouses <- create_spouse_file(agta_genealogies, AGTA_BAYAKA, agta_coef) 

# Save datafiles
write.csv(bayaka_spouses, "/Users/inezd/Documents/Science/Raute/Chapter_3/R_Files/bayaka_spouses.csv")
write.csv(agta_spouses, "/Users/inezd/Documents/Science/Raute/Chapter_3/R_Files/agta_spouses.csv")


##################################### ELIGIBLE POPULATION #####################################

### For the Agta

# Add age
agta_data <- agta_data %>%
  mutate(age = 2014 - `Birth Year`)

# Create file
agta_eligible <- as.data.frame(t(combn(agta_data$ID, 2))) %>% # make combinations of all individuals
  rename(X1 = V1, X2 = V2) %>% # rename
  left_join(agta_data %>% select(ID, Sex, age), by = c("X1" = "ID")) %>% # add age and sex of one individual
  left_join(agta_data %>% select(ID, Sex, age), by = c("X2" = "ID"), suffix = c(".x", ".y")) %>% # add age and sex of second individual
  filter(Sex.x != Sex.y, !is.na(age.x), !is.na(age.y), age.x > 16, age.y > 16) %>% # select only opposite sex pairs, have to have age, and older than 16
  left_join(AGTA_BAYAKA %>% select(X1 = old_FID, new_IID), by = "X1") %>% # add new id
  left_join(AGTA_BAYAKA %>% select(X2 = old_FID, new_IID), by = "X2", suffix = c(".x", ".y")) %>% # add new id
  mutate(pair = paste(pmin(new_IID.x, new_IID.y), pmax(new_IID.x, new_IID.y), sep = "-")) %>%
  left_join(agta_coef, by = "pair") %>%
  select(pair, kinship) %>%
  distinct() %>%
  filter(kinship < 0.1768) #2422, mean is 0.02020264


### For the BaYaka
bayaka_eligible <- as.data.frame(t(combn(bayaka_data$short_ID[!bayaka_data$short_ID == 'X236'], 2))) %>% # make combinations of all individuals
  rename(X1 = V1, X2 = V2) %>% # rename
  distinct() %>%
  left_join(bayaka_data %>% select(X1 = short_ID, Gender, LifeStage), by = "X1") %>% # add age and sex of one individual
  left_join(bayaka_data %>% select(X2 = short_ID, Gender, LifeStage), by = "X2", suffix = c(".x", ".y")) %>% # add age and sex of second individual
  filter(Gender.x != Gender.y, grepl("adult", LifeStage.x), grepl("adult", LifeStage.y)) %>%
  left_join(AGTA_BAYAKA %>% select(X1 = old_FID, new_IID), by = "X1") %>% # add new id
  left_join(AGTA_BAYAKA %>% select(X2 = old_FID, new_IID), by = "X2", suffix = c(".x", ".y")) %>% # add new id
  mutate(pair = paste(pmin(new_IID.x, new_IID.y), pmax(new_IID.x, new_IID.y), sep = "-")) %>%
  left_join(bayaka_coef, by = "pair") %>%
  select(pair, kinship) %>%
  distinct() %>%
  filter(kinship < 0.1768) #49, mean is 0.02500855



##################################### TOTAL POPULATION #####################################

### For the Agta
agta_total <- agta_coef %>% # relatedness file
  left_join(agta_spouses %>% select(pair, type), by = "pair") %>% # add spouses label to the datafile
  mutate(type = ifelse(is.na(type), "non-spouse", type)) %>% # add the label 'non-spouse' to those who are not spouses
  select(pair, kinship, type) %>%
  distinct() # mean = 0.02539593, n = 9755

### For the BaYaka
bayaka_total <- bayaka_coef %>% # relatedness file
  left_join(bayaka_spouses %>% select(pair, type), by = "pair") %>% # add spouses label to the datafile
  mutate(type = ifelse(is.na(type), "non-spouse", type)) %>% # add the label 'non-spouse' to those who are not spouses
  select(pair, kinship, type) %>%
  distinct() # mean = 0.04604954, n = 166


##################################### PERMUTE RELATEDNESS #####################################

# Use function from script 'Function_Permute_relatedness.R'

bayaka_permuted <- permute_relatedness(bayaka_total)
agta_permuted <- permute_relatedness(agta_total)
































