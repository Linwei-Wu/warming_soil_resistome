library(tidyr)
# Genomic Proximity Analysis in R

# this script tests the hypothesis that ARGs are genomically linked
#              to stress response genes (High Temp Tol, N Assim, C Degrad).

# --- 1. Load Libraries ---
setwd("/home/linwei/nws0920/22_genomic_proximity")

library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(broom) # For tidying model outputs

# --- 2. Configuration ---
# Define the gene sets of interest
TARGET_GENE_SETS <- c("Tol_high_temp", "N_assimilation", "complex_C_degrad", "MGE")
# Define the window size (number of genes upstream and downstream of the ARG)
WINDOW_SIZE <- 10
# Number of permutations for the null model
N_PERMUTATIONS <- 1000

# --- 3. Load and Prepare Data ---
cat("Loading gene annotation data...\n")

gene_df <- read_csv("gene_df_strains.csv")
gene_df <- gene_df %>% filter(phy_group == "Actinobacteria") 
head(gene_df) #see how the table is structured below

#> head(gene_df) 
## A tibble: 6 × 9
#gene_id        start   end strand info  genome_id contig_id gene_set phy_group
#<chr>          <dbl> <dbl>  <dbl> <chr> <chr>     <chr>     <chr>    <chr>    
#  1 Q318.scaffold…   213  1820     -1 ID=1… Q318      Q318.sca… simple_… Actinoba…
#2 Q318.scaffold…  1907  2401      1 ID=1… Q318      Q318.sca… houseke… Actinoba…
#3 Q318.scaffold…  2560  3414     -1 ID=1… Q318      Q318.sca… others   Actinoba…
#4 Q318.scaffold…  3558  3740     -1 ID=1… Q318      Q318.sca… others   Actinoba…
#5 Q318.scaffold…  3868  4455     -1 ID=1… Q318      Q318.sca… others   Actinoba…
#6 Q318.scaffold…  4452  5267     -1 ID=1… Q318      Q318.sca… others   Actinoba…




# Data cleaning and preparation
# Arrange genes by their physical order on the contig
cat("Preparing data by sorting genes...\n")
gene_df_sorted <- gene_df %>%
  arrange(genome_id, contig_id, start) %>%
  group_by(genome_id, contig_id) %>%
  mutate(gene_index = row_number()) # Add a gene index within each contig

# --- 4. Function for Proximity Analysis ---
# This function finds neighbors for a focal gene set (e.g., ARGs)
find_neighbors <- function(df, focal_set = "ARG", target_sets, window = 10) {
  
  focal_indices <- which(df$gene_set == focal_set)
  
  if (length(focal_indices) == 0) {
    return(tibble()) # Return empty tibble if no focal genes
  }
  
  neighbor_list <- list()
  
  for (i in focal_indices) {
    # Define the window boundaries
    start_window <- max(1, i - window)
    end_window <- min(nrow(df), i + window)
    
    window_df <- df[start_window:end_window, ]
    
    # Exclude the focal gene itself from its neighborhood
    window_df <- window_df %>% filter(gene_index != df$gene_index[i])
    
    # Find which target sets are present in the window
    present_targets <- intersect(target_sets, unique(window_df$gene_set))
    
    if (length(present_targets) > 0) {
      neighbor_list[[as.character(i)]] <- tibble(
        focal_gene_id = df$gene_id[i],
        neighbor_gene_set = present_targets
      )
    }
  }
  
  if (length(neighbor_list) > 0) {
    return(bind_rows(neighbor_list))
  } else {
    return(tibble())
  }
}

# --- 5. Observed Co-occurrence Analysis ---
cat("Calculating observed co-occurrence of ARGs and target gene sets...\n")
# We apply the function to each contig separately
observed_results <- gene_df_sorted %>%
  group_by(genome_id, contig_id) %>%
  nest() %>%
  mutate(
    neighbors = map(data, ~find_neighbors(., focal_set = "ARG", target_sets = TARGET_GENE_SETS, window = WINDOW_SIZE))
  ) %>%
  unnest(neighbors) %>%
  select(genome_id, contig_id, focal_gene_id, neighbor_gene_set)

phygroup <- gene_df$phy_group[match(observed_results$genome_id, gene_df$genome_id)]
observed_results <- data.frame(observed_results, phy_group = phygroup)
write_csv(observed_results, "observed_co-occur_strains_actino.csv")

# Summarize the observed counts
observed_counts <- observed_results %>%
  group_by(neighbor_gene_set) %>%
  summarise(observed_count = n(), .groups = 'drop')

cat("Observed co-occurrence counts:\n")
print(observed_counts)

# --- 6. Permutation Test (Null Model) ---
cat(paste0("Running permutation test with ", N_PERMUTATIONS, " iterations...\n"))

# This function shuffles the 'gene_set' column within a contig
permute_and_analyze <- function(df_nested, .id) {
  df_nested %>%
    mutate(
      shuffled_data = map(data, ~ .x %>% mutate(gene_set = sample(gene_set))),
      shuffled_neighbors = map(shuffled_data, ~find_neighbors(., focal_set = "ARG", target_sets = TARGET_GENE_SETS, window = WINDOW_SIZE))
    ) %>%
    unnest(shuffled_neighbors) %>%
    group_by(neighbor_gene_set) %>%
    summarise(permuted_count = n(), .groups = 'drop')
}


# We use map from purrr to iterate. Replace with future_map for parallel execution.
permuted_results <- map_dfr(
  1:N_PERMUTATIONS,
  ~ permute_and_analyze(gene_df_sorted %>% group_by(genome_id, contig_id) %>% nest(), .id = .x),
  .id = "permutation_id"
)

# Summarize permutation results
permuted_summary <- permuted_results %>%
  group_by(neighbor_gene_set) %>%
  summarise(
    mean_perm_count = mean(permuted_count),
    sd_perm_count = sd(permuted_count),
    .groups = 'drop'
  )


# --- 7. Calculate Enrichment and P-values ---
cat("Calculating enrichment and p-values...\n")
final_results <- observed_counts %>%
  left_join(permuted_summary, by = "neighbor_gene_set") %>%
  mutate(
    # Fold enrichment = Observed / Mean of Random
    fold_enrichment = observed_count / mean_perm_count,
    # Calculate empirical p-value
    p_value = map_dbl(neighbor_gene_set, ~ {
      obs_c <- observed_count[neighbor_gene_set == .x]
      perm_counts <- permuted_results %>% filter(neighbor_gene_set == .x) %>% pull(permuted_count)
      # P-value is the proportion of permuted counts >= observed count
      (sum(perm_counts >= obs_c) + 1) / (N_PERMUTATIONS + 1)
    })
  ) %>% 
  arrange(p_value)

cat("Final enrichment results:\n")
print(final_results)

write_csv(final_results, "final_results_strains_actino.csv")



