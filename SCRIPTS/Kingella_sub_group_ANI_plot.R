
master_ANI_pangenome_ID_order <- as.data.frame(levels(long_df$key))

# Subset keys dataframe
Kingela_sub_group_keys <- master_ANI_pangenome_ID_order[124:150,, drop = FALSE]


# Convert the dataframe column to a vector
Kingela_sub_group_pangenomeIDs <- Kingela_sub_group_keys[,1]

# Subset your matrix
Kingela_sub_group_long_df <- long_df %>% 
  filter(key %in% Kingela_sub_group_pangenomeIDs) %>% 
  filter(y %in% Kingela_sub_group_pangenomeIDs) %>% 
  droplevels()


Kingela_sub_group_order <- long_df %>% 
  filter(key %in% Kingela_sub_group_pangenomeIDs) %>% 
  filter(y %in% Kingela_sub_group_pangenomeIDs) %>% 
  droplevels() %>% 
  select(key)

Kingela_sub_group_order_list <- levels(Kingela_sub_group_order$key)


colors_Kingela_sub_group <- as.data.frame(levels(Kingela_sub_group_long_df$key))
colors_Kingela_sub_group <- colors_Kingela_sub_group %>% 
  rename(pangenomeIDs = "levels(Kingela_sub_group_long_df$key)") %>% 
  mutate(colors = case_when(grepl('K_negevensis', pangenomeIDs) ~ 'blue',
                            grepl('N_arctica', pangenomeIDs) ~ 'black',
                            grepl('K_denitrificans', pangenomeIDs) ~ 'magenta',
                            grepl('K_oralis', pangenomeIDs) ~ 'darkgreen',
                            grepl('N_shayeganii', pangenomeIDs) ~ 'darkorange',
                            grepl('N_sp', pangenomeIDs) ~ 'red',
                            grepl('K_bonacorsii', pangenomeIDs) ~ 'palegreen',
                            grepl('N_montereyensis', pangenomeIDs) ~ 'black'))


# Deal with duplicate strains K_negevensis_str_NA
#K_negevensis_str_NA_id_GCA_000751855_1 K_negevensis_str_NA_A_id_GCA_000751855_1
#K_negevensis_str_NA_id_GCA_900182485_2 K_negevensis_str_NA_B_id_GCA_900182485_2
colors_Kingela_sub_group$pangenomeIDs <- gsub("K_negevensis_str_NA_id_GCA_000751855_1", "K_negevensis_str_NA_A_id_GCA_000751855_1", colors_Kingela_sub_group$pangenomeIDs)
colors_Kingela_sub_group$pangenomeIDs <- gsub("K_negevensis_str_NA_id_GCA_900182485_2", "K_negevensis_str_NA_B_id_GCA_900182485_2", colors_Kingela_sub_group$pangenomeIDs)
Kingela_sub_group_long_df$key <- gsub("K_negevensis_str_NA_id_GCA_000751855_1", "K_negevensis_str_NA_A_id_GCA_000751855_1", Kingela_sub_group_long_df$key)
Kingela_sub_group_long_df$y <- gsub("K_negevensis_str_NA_id_GCA_000751855_1", "K_negevensis_str_NA_A_id_GCA_000751855_1", Kingela_sub_group_long_df$y)
Kingela_sub_group_long_df$key <- gsub("K_negevensis_str_NA_id_GCA_900182485_2", "K_negevensis_str_NA_B_id_GCA_900182485_2", Kingela_sub_group_long_df$key)
Kingela_sub_group_long_df$y <- gsub("K_negevensis_str_NA_id_GCA_900182485_2", "K_negevensis_str_NA_B_id_GCA_900182485_2", Kingela_sub_group_long_df$y)
Kingela_sub_group_order_list <- gsub("K_negevensis_str_NA_id_GCA_000751855_1", "K_negevensis_str_NA_A_id_GCA_000751855_1", Kingela_sub_group_order_list)
Kingela_sub_group_order_list <- gsub("K_negevensis_str_NA_id_GCA_900182485_2", "K_negevensis_str_NA_B_id_GCA_900182485_2", Kingela_sub_group_order_list)

#K_negevensis_str_SW7208426
#K_negevensis_str_SW7208426_id_GCA_900177895_1 K_negevensis_str_SW7208426_A_id_GCA_900177895_1
#K_negevensis_str_SW7208426_id_GCA_030177895_1 K_negevensis_str_SW7208426_B_id_GCA_030177895_1
colors_Kingela_sub_group$pangenomeIDs <- gsub("K_negevensis_str_SW7208426_id_GCA_900177895_1", "K_negevensis_str_SW7208426_A_id_GCA_900177895_1", colors_Kingela_sub_group$pangenomeIDs)
colors_Kingela_sub_group$pangenomeIDs <- gsub("K_negevensis_str_SW7208426_id_GCA_030177895_1", "K_negevensis_str_SW7208426_B_id_GCA_030177895_1", colors_Kingela_sub_group$pangenomeIDs)

Kingela_sub_group_long_df$key <- gsub("K_negevensis_str_SW7208426_id_GCA_900177895_1", "K_negevensis_str_SW7208426_A_id_GCA_900177895_1", Kingela_sub_group_long_df$key)
Kingela_sub_group_long_df$y <- gsub("K_negevensis_str_SW7208426_id_GCA_900177895_1", "K_negevensis_str_SW7208426_A_id_GCA_900177895_1", Kingela_sub_group_long_df$y)

Kingela_sub_group_long_df$key <- gsub("K_negevensis_str_SW7208426_id_GCA_030177895_1", "K_negevensis_str_SW7208426_B_id_GCA_030177895_1", Kingela_sub_group_long_df$key)
Kingela_sub_group_long_df$y <- gsub("K_negevensis_str_SW7208426_id_GCA_030177895_1", "K_negevensis_str_SW7208426_B_id_GCA_030177895_1", Kingela_sub_group_long_df$y)

Kingela_sub_group_order_list <- gsub("K_negevensis_str_SW7208426_id_GCA_900177895_1", "K_negevensis_str_SW7208426_A_id_GCA_900177895_1", Kingela_sub_group_order_list)
Kingela_sub_group_order_list <- gsub("K_negevensis_str_SW7208426_id_GCA_030177895_1", "K_negevensis_str_SW7208426_B_id_GCA_030177895_1", Kingela_sub_group_order_list)



#N_arctica_str_KH1503
#N_arctica_str_KH1503_id_GCA_022870905_1 N_arctica_str_KH1503_A_id_GCA_022870905_1
#N_arctica_str_KH1503_id_GCA_001027865_1 N_arctica_str_KH1503_B_id_GCA_001027865_1
colors_Kingela_sub_group$pangenomeIDs <- gsub("N_arctica_str_KH1503_id_GCA_022870905_1", "N_arctica_str_KH1503_A_id_GCA_022870905_1", colors_Kingela_sub_group$pangenomeIDs)
colors_Kingela_sub_group$pangenomeIDs <- gsub("N_arctica_str_KH1503_id_GCA_001027865_1", "N_arctica_str_KH1503_B_id_GCA_001027865_1", colors_Kingela_sub_group$pangenomeIDs)
Kingela_sub_group_long_df$key <- gsub("N_arctica_str_KH1503_id_GCA_022870905_1", "N_arctica_str_KH1503_A_id_GCA_022870905_1", Kingela_sub_group_long_df$key)
Kingela_sub_group_long_df$y <- gsub("N_arctica_str_KH1503_id_GCA_022870905_1", "N_arctica_str_KH1503_A_id_GCA_022870905_1", Kingela_sub_group_long_df$y)
Kingela_sub_group_long_df$key <- gsub("N_arctica_str_KH1503_id_GCA_001027865_1", "N_arctica_str_KH1503_B_id_GCA_001027865_1", Kingela_sub_group_long_df$key)
Kingela_sub_group_long_df$y <- gsub("N_arctica_str_KH1503_id_GCA_001027865_1", "N_arctica_str_KH1503_B_id_GCA_001027865_1", Kingela_sub_group_long_df$y)
Kingela_sub_group_order_list <- gsub("N_arctica_str_KH1503_id_GCA_022870905_1", "N_arctica_str_KH1503_A_id_GCA_022870905_1", Kingela_sub_group_order_list)
Kingela_sub_group_order_list <- gsub("N_arctica_str_KH1503_id_GCA_001027865_1", "N_arctica_str_KH1503_B_id_GCA_001027865_1", Kingela_sub_group_order_list)



# remove GCA id from names
colors_Kingela_sub_group$pangenomeIDs <- gsub("_id_.*", "", colors_Kingela_sub_group$pangenomeIDs)
Kingela_sub_group_long_df$key <- gsub("_id_.*", "", Kingela_sub_group_long_df$key)
Kingela_sub_group_long_df$y <- gsub("_id_.*", "", Kingela_sub_group_long_df$y)
Kingela_sub_group_order_list <- gsub("_id_.*", "", Kingela_sub_group_order_list)


# Reorder levels of factor key and y
Kingela_sub_group_long_df$key <- as.factor(Kingela_sub_group_long_df$key)
Kingela_sub_group_long_df$y <- as.factor(Kingela_sub_group_long_df$y)
Kingela_sub_group_long_df$key <- factor(Kingela_sub_group_long_df$key, levels = Kingela_sub_group_order_list)  
Kingela_sub_group_long_df$y <- factor(Kingela_sub_group_long_df$y, levels = Kingela_sub_group_order_list)  


# rename z to ANI
Kingela_sub_group_long_df <- Kingela_sub_group_long_df %>% 
  rename(ANI = "z")

# Then, create the plot:
plot <- ggplot(Kingela_sub_group_long_df, aes(key, y)) +
  geom_tile(aes(fill = ANI)) + 
  scale_fill_gradient2(low = "darkblue",
                       mid = "white",
                       high = "darkred",
                       midpoint = 0.95,
                       limits =c(0.90, 1),
                       na.value="darkblue") +
  ylab("") +  
  xlab("") + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1, color = colors_Kingela_sub_group$colors),
        axis.text.y = element_text(size = 10, color = colors_Kingela_sub_group$colors),
        axis.ticks=element_blank())

ggsave(file = "/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/ANI/G_0636/G_0636_ANI_heatmap_colors_Kingela_sub_group_test.pdf", plot = plot, width = 11, height = 8, limitsize = FALSE)

