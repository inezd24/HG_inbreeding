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

##################################### TRANSFORM DATA #####################################

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
    mutate(type = 'spouse') # add variable that shows these are spouses
  return(out_file) # return output file
}
bayaka_spouses <- create_spouse_file(bayaka_data, AGTA_BAYAKA, bayaka_coef)
agta_spouses <- create_spouse_file(agta_genealogies, AGTA_BAYAKA, agta_coef) 

# Those datafiles are not finished yet, first we have to make an 'individuals' file
make_individual_file <- function(spouses_file){
  ind_file <- spouses_file %>%
    select(kinship, sample_1, sample_2) %>%
    pivot_longer(!kinship, values_to = 'sample_1') %>%
    distinct(sample_1)
  return(ind_file)
}
bayaka_spouse_ind <- make_individual_file(bayaka_spouses)
agta_spouse_ind <- make_individual_file(agta_spouses)

# Finish previous data files
finish_spouses_file <- function(spouses_file){
  finished_file <- spouses_file %>%
    select(pair, kinship, type) %>%
    distinct() 
  return(finished_file)
}
bayaka_spouses <- finish_spouses_file(bayaka_spouses)
agta_spouses <- finish_spouses_file(agta_spouses)

# Save datafiles
write.csv(bayaka_spouses, "/Users/inezd/Documents/Science/Raute/Chapter_3/R_Files/bayaka_spouses.csv")
write.csv(agta_spouses, "/Users/inezd/Documents/Science/Raute/Chapter_3/R_Files/agta_spouses.csv")



















