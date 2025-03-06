# Script written by Inez Derkx, edited March 2025
# This script is a function that allows you to permute relatedness in the datafile structure used in the paper '', Derkx et al. 2025. 
# This function requires a datafile that has the following variables and assumes you want to run 1000 permutations with a seed of 123:
      # kinship: this is the dyadic relatedness between a pair of individuals. It can be any kinship measure you desire, as long as it's numerical. 
      # type: this is the type of relationship of the dyad. This variable is a character/string. 
      # It is assumed that the datafile contains pairs of individuals who are spouses (named 'spouse') and those who aren't spouses (named 'non-spouse'). 

permute_relatedness <- function(datafile, num_permutations = 1000, seed = 123) {
  
  # For reproducibility
  set.seed(seed)  
  
  # Save the observed mean among spouses
  observed_mean <- mean(datafile$sum[datafile$type == 'spouse'])
  
  # Set a number of permutations
  permuted_means <- replicate(num_permutations, {
    
    # Shuffle the kinship coefficients
    shuffled_relationship <- sample(datafile$kinship)
    
    # After shuffling, recalculate the mean kinship coefficient for those named spouse
    mean(shuffled_relationship[datafile$type == 'spouse'])
  })
  
  # Save in new datafile
  permuted_means_df <- data.frame(sim = permuted_means)
  
  # Calculate statistics
  mean_distribution <- mean(permuted_means_df$sim)
  sd_distribution <- sd(permuted_means_df$sim)
  z_score <- (observed_mean - mean_distribution) / sd_distribution
  p_value <- pnorm(z_score)
  
  # Return all as list
  return(list(
    observed_mean = observed_mean,
    mean_distribution = mean_distribution,
    sd_distribution = sd_distribution,
    z_score = z_score,
    p_value = p_value,
    permuted_means_df = permuted_means_df
  ))
}
