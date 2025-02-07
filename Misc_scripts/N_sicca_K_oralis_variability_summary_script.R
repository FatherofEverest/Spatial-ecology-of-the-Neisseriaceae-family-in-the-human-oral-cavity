#!/usr/bin/env Rscript

library("dplyr")

##### N.sicca ##### 

N_sicca_df <- read.table("/workspace/jmarkwelchlab/P_0003_Neisseriaceae/16_VARIABILITY/N_sicca_str_VK64_id_GCA_000260655_1-BM-SNVs.txt", header = TRUE, stringsAsFactors = TRUE, sep = "\t")

N_sicca_PP_df <- read.table("/workspace/jmarkwelchlab/P_0003_Neisseriaceae/16_VARIABILITY/N_sicca_str_VK64_id_GCA_000260655_1-PP-SNVs.txt", header = TRUE, stringsAsFactors = TRUE, sep = "\t")

N_sicca_genome_length_kb <- 2638.449  # Replace with actual genome length in kb

N_sicca_BM_summary_metrics <- N_sicca_df %>%
  filter(departure_from_reference > 0) %>% 
  filter(coverage >= 10) %>%
  group_by(sample_id) %>%
  summarise(
    mean_depart_ref = mean(departure_from_reference, na.rm = TRUE),
    mean_depart_con = mean(departure_from_consensus, na.rm = TRUE),
    mean_kl_divergence_raw = mean(kullback_leibler_divergence_raw, na.rm = TRUE),
    mean_entropy = mean(entropy, na.rm = TRUE),
    mean_kl_divergence_normalized = mean(kullback_leibler_divergence_normalized, na.rm = TRUE),
    weighted_average_entropy = sum(entropy * coverage, na.rm = TRUE) / sum(coverage, na.rm = TRUE),
    snv_density = n_distinct(unique_pos_identifier) / N_sicca_genome_length_kb)
N_sicca_BM_summary_metrics <- N_sicca_BM_summary_metrics %>% 
  mutate(species = as.factor("N_sicca")) %>% 
  mutate(site = as.factor(case_when(grepl('PP', sample_id) ~ 'SUPP',
                          grepl('BM', sample_id) ~ 'BM')))

N_sicca_PP_summary_metrics <- N_sicca_PP_df %>%
  filter(departure_from_reference > 0) %>% 
  filter(coverage >= 10) %>%
  group_by(sample_id) %>%
  summarise(
    mean_depart_ref = mean(departure_from_reference, na.rm = TRUE),
    mean_depart_con = mean(departure_from_consensus, na.rm = TRUE),
    mean_kl_divergence_raw = mean(kullback_leibler_divergence_raw, na.rm = TRUE),
    mean_entropy = mean(entropy, na.rm = TRUE),
    mean_kl_divergence_normalized = mean(kullback_leibler_divergence_normalized, na.rm = TRUE),
    weighted_average_entropy = sum(entropy * coverage, na.rm = TRUE) / sum(coverage, na.rm = TRUE),
    snv_density = n_distinct(unique_pos_identifier) / N_sicca_genome_length_kb)
N_sicca_PP_summary_metrics <- N_sicca_PP_summary_metrics %>% 
  mutate(species = as.factor("N_sicca")) %>% 
  mutate(site = as.factor(case_when(grepl('PP', sample_id) ~ 'SUPP',
                          grepl('BM', sample_id) ~ 'BM')))

N_sicca_combined_metrics <- bind_rows(N_sicca_BM_summary_metrics, N_sicca_PP_summary_metrics)


#####  K.oralis ##### 

K_oralis_PP_df <- read.table("/workspace/jmarkwelchlab/P_0003_Neisseriaceae/16_VARIABILITY/K_oralis_str_ATCC_51147_id_GCA_000160435_1-PP-SNVs.txt", header = TRUE, stringsAsFactors = TRUE, sep = "\t")

K_oralis_BM_df <- read.table("/workspace/jmarkwelchlab/P_0003_Neisseriaceae/16_VARIABILITY/K_oralis_str_ATCC_51147_id_GCA_000160435_1-BM-SNVs.txt", header = TRUE, stringsAsFactors = TRUE, sep = "\t")

K_oralis_genome_length_kb <- 2406.675  # Replace with actual genome length in kb


K_oralis_PP_summary_metrics <- K_oralis_PP_df %>%
  filter(departure_from_reference > 0) %>% 
  filter(coverage >= 10) %>%
  group_by(sample_id) %>%
  summarise(
    mean_depart_ref = mean(departure_from_reference, na.rm = TRUE),
    mean_depart_con = mean(departure_from_consensus, na.rm = TRUE),
    mean_kl_divergence_raw = mean(kullback_leibler_divergence_raw, na.rm = TRUE),
    mean_entropy = mean(entropy, na.rm = TRUE),
    mean_kl_divergence_normalized = mean(kullback_leibler_divergence_normalized, na.rm = TRUE),
    weighted_average_entropy = sum(entropy * coverage, na.rm = TRUE) / sum(coverage, na.rm = TRUE),
    snv_density = n_distinct(unique_pos_identifier) / K_oralis_genome_length_kb)
K_oralis_PP_summary_metrics <- K_oralis_PP_summary_metrics %>% 
  mutate(species = as.factor("K_oralis")) %>% 
  mutate(site = as.factor(case_when(grepl('PP', sample_id) ~ 'SUPP',
                          grepl('BM', sample_id) ~ 'BM')))


K_oralis_BM_summary_metrics <- K_oralis_BM_df %>%
  filter(departure_from_reference > 0) %>% 
  filter(coverage >= 10) %>%
  group_by(sample_id) %>%
  summarise(
    mean_depart_ref = mean(departure_from_reference, na.rm = TRUE),
    mean_depart_con = mean(departure_from_consensus, na.rm = TRUE),
    mean_kl_divergence_raw = mean(kullback_leibler_divergence_raw, na.rm = TRUE),
    mean_entropy = mean(entropy, na.rm = TRUE),
    mean_kl_divergence_normalized = mean(kullback_leibler_divergence_normalized, na.rm = TRUE),
    weighted_average_entropy = sum(entropy * coverage, na.rm = TRUE) / sum(coverage, na.rm = TRUE),
    snv_density = n_distinct(unique_pos_identifier) / K_oralis_genome_length_kb)
K_oralis_BM_summary_metrics <- K_oralis_BM_summary_metrics %>% 
  mutate(species = as.factor("K_oralis")) %>% 
  mutate(site = as.factor(case_when(grepl('PP', sample_id) ~ 'SUPP',
                          grepl('BM', sample_id) ~ 'BM')))

K_oralis_combined_metrics <- bind_rows(K_oralis_PP_summary_metrics, K_oralis_BM_summary_metrics)


final_data <- bind_rows(K_oralis_combined_metrics,N_sicca_combined_metrics)

write.table(final_data, "/workspace/jmarkwelchlab/P_0003_Neisseriaceae/16_VARIABILITY/N_sicca_and_K_oralis_variability_summary.txt", sep = "\t", quote = FALSE, row.names = FALSE)

