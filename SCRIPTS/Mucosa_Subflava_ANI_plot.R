
master_ANI_pangenome_ID_order <- as.data.frame(levels(long_df$key))

# Subset keys dataframe
Mucosa_Subflava_keys <- master_ANI_pangenome_ID_order[388:537,, drop = FALSE]


# Convert the dataframe column to a vector
Mucosa_Subflava_pangenomeIDs <- Mucosa_Subflava_keys[,1]

# Subset your matrix
Mucosa_Subflava_long_df <- long_df %>% 
  filter(key %in% Mucosa_Subflava_pangenomeIDs) %>% 
  filter(y %in% Mucosa_Subflava_pangenomeIDs) %>% 
  droplevels()


Mucosa_Subflava_order <- long_df %>% 
  filter(key %in% Mucosa_Subflava_pangenomeIDs) %>% 
  filter(y %in% Mucosa_Subflava_pangenomeIDs) %>% 
  droplevels() %>% 
  select(key)

Mucosa_Subflava_order_list <- levels(Mucosa_Subflava_order$key)


colors_Mucosa_Subflava <- as.data.frame(levels(Mucosa_Subflava_long_df$key))
colors_Mucosa_Subflava <- colors_Mucosa_Subflava %>% 
  rename(pangenomeIDs = "levels(Mucosa_Subflava_long_df$key)") %>% 
  mutate(colors = case_when(grepl('N_sp', pangenomeIDs) ~ 'blue',
                            grepl('N_cinerea', pangenomeIDs) ~ 'cyan',
                            grepl('N_lactamica', pangenomeIDs) ~ 'red',
                            grepl('N_polysaccharea', pangenomeIDs) ~ 'steelblue',
                            grepl('N_flavescens', pangenomeIDs) ~ 'darkgreen',
                            grepl('N_elongata', pangenomeIDs) ~ 'green',
                            grepl('N_subflava', pangenomeIDs) ~ 'magenta',
                            grepl('N_oralis', pangenomeIDs) ~ 'olivedrab',
                            grepl('N_perflava', pangenomeIDs) ~ 'goldenrod1',
                            grepl('N_sicca', pangenomeIDs) ~ 'purple',
                            grepl('N_mucosa', pangenomeIDs) ~ 'darkred',
                            grepl('N_macacae', pangenomeIDs) ~ 'black',))

# remove GCA id from names
#N_sp_str_F0314 dulpicate
#N_sp_str_F0314_id_GCA_000090875_1 N_sp_str_F0314_A_id_GCA_000090875_1
#N_sp_str_F0314_id_GCA_005886145_1 N_sp_str_F0314_B_id_GCA_005886145_1
colors_Mucosa_Subflava$pangenomeIDs <- gsub("N_sp_str_F0314_id_GCA_000090875_1", "N_sp_str_F0314_A_id_GCA_000090875_1", colors_Mucosa_Subflava$pangenomeIDs)
colors_Mucosa_Subflava$pangenomeIDs <- gsub("N_sp_str_F0314_id_GCA_005886145_1", "N_sp_str_F0314_B_id_GCA_005886145_1", colors_Mucosa_Subflava$pangenomeIDs)

Mucosa_Subflava_long_df$key <- gsub("N_sp_str_F0314_id_GCA_000090875_1", "N_sp_str_F0314_A_id_GCA_000090875_1", Mucosa_Subflava_long_df$key)
Mucosa_Subflava_long_df$y <- gsub("N_sp_str_F0314_id_GCA_000090875_1", "N_sp_str_F0314_A_id_GCA_000090875_1", Mucosa_Subflava_long_df$y)

Mucosa_Subflava_long_df$key <- gsub("N_sp_str_F0314_id_GCA_005886145_1", "N_sp_str_F0314_B_id_GCA_005886145_1", Mucosa_Subflava_long_df$key)
Mucosa_Subflava_long_df$y <- gsub("N_sp_str_F0314_id_GCA_005886145_1", "N_sp_str_F0314_B_id_GCA_005886145_1", Mucosa_Subflava_long_df$y)

Mucosa_Subflava_order_list <- gsub("N_sp_str_F0314_id_GCA_000090875_1", "N_sp_str_F0314_A_id_GCA_000090875_1", Mucosa_Subflava_order_list)
Mucosa_Subflava_order_list <- gsub("N_sp_str_F0314_id_GCA_005886145_1", "N_sp_str_F0314_B_id_GCA_005886145_1", Mucosa_Subflava_order_list)


#N_sicca_str_ATCC_29256 dulpicate _id_GCA_000174655_1 _id_GCA_019334765_1
colors_Mucosa_Subflava$pangenomeIDs <- gsub("N_sicca_str_ATCC_29256_id_GCA_000174655_1", "N_sicca_str_ATCC_29256_A_id_GCA_000174655_1", colors_Mucosa_Subflava$pangenomeIDs)
colors_Mucosa_Subflava$pangenomeIDs <- gsub("N_sicca_str_ATCC_29256_id_GCA_019334765_1", "N_sicca_str_ATCC_29256_B_id_GCA_019334765_1", colors_Mucosa_Subflava$pangenomeIDs)

Mucosa_Subflava_long_df$key <- gsub("N_sicca_str_ATCC_29256_id_GCA_000174655_1", "N_sicca_str_ATCC_29256_A_id_GCA_000174655_1", Mucosa_Subflava_long_df$key)
Mucosa_Subflava_long_df$y <- gsub("N_sicca_str_ATCC_29256_id_GCA_000174655_1", "N_sicca_str_ATCC_29256_A_id_GCA_000174655_1", Mucosa_Subflava_long_df$y)

Mucosa_Subflava_long_df$key <- gsub("N_sicca_str_ATCC_29256_id_GCA_019334765_1", "N_sicca_str_ATCC_29256_B_id_GCA_019334765_1", Mucosa_Subflava_long_df$key)
Mucosa_Subflava_long_df$y <- gsub("N_sicca_str_ATCC_29256_id_GCA_019334765_1", "N_sicca_str_ATCC_29256_B_id_GCA_019334765_1", Mucosa_Subflava_long_df$y)

Mucosa_Subflava_order_list <- gsub("N_sicca_str_ATCC_29256_id_GCA_000174655_1", "N_sicca_str_ATCC_29256_A_id_GCA_000174655_1", Mucosa_Subflava_order_list)
Mucosa_Subflava_order_list <- gsub("N_sicca_str_ATCC_29256_id_GCA_019334765_1", "N_sicca_str_ATCC_29256_B_id_GCA_019334765_1", Mucosa_Subflava_order_list)

#N_macacae_str_ATCC_33926 _id_GCA_022749495_1 _id_GCA_000220865_1
colors_Mucosa_Subflava$pangenomeIDs <- gsub("N_macacae_str_ATCC_33926_id_GCA_022749495_1", "N_macacae_str_ATCC_33926_A_id_GCA_022749495_1", colors_Mucosa_Subflava$pangenomeIDs)
colors_Mucosa_Subflava$pangenomeIDs <- gsub("N_macacae_str_ATCC_33926_id_GCA_000220865_1", "N_macacae_str_ATCC_33926_B_id_GCA_000220865_1", colors_Mucosa_Subflava$pangenomeIDs)

Mucosa_Subflava_long_df$key <- gsub("N_macacae_str_ATCC_33926_id_GCA_022749495_1", "N_macacae_str_ATCC_33926_A_id_GCA_022749495_1", Mucosa_Subflava_long_df$key)
Mucosa_Subflava_long_df$y <- gsub("N_macacae_str_ATCC_33926_id_GCA_022749495_1", "N_macacae_str_ATCC_33926_A_id_GCA_022749495_1", Mucosa_Subflava_long_df$y)

Mucosa_Subflava_long_df$key <- gsub("N_macacae_str_ATCC_33926_id_GCA_000220865_1", "N_macacae_str_ATCC_33926_B_id_GCA_000220865_1", Mucosa_Subflava_long_df$key)
Mucosa_Subflava_long_df$y <- gsub("N_macacae_str_ATCC_33926_id_GCA_000220865_1", "N_macacae_str_ATCC_33926_B_id_GCA_000220865_1", Mucosa_Subflava_long_df$y)

Mucosa_Subflava_order_list <- gsub("N_macacae_str_ATCC_33926_id_GCA_022749495_1", "N_macacae_str_ATCC_33926_A_id_GCA_022749495_1", Mucosa_Subflava_order_list)
Mucosa_Subflava_order_list <- gsub("N_macacae_str_ATCC_33926_id_GCA_000220865_1", "N_macacae_str_ATCC_33926_B_id_GCA_000220865_1", Mucosa_Subflava_order_list)


colors_Mucosa_Subflava$pangenomeIDs <- gsub("_id_.*", "", colors_Mucosa_Subflava$pangenomeIDs)
Mucosa_Subflava_long_df$key <- gsub("_id_.*", "", Mucosa_Subflava_long_df$key)
Mucosa_Subflava_long_df$y <- gsub("_id_.*", "", Mucosa_Subflava_long_df$y)
Mucosa_Subflava_order_list <- gsub("_id_.*", "", Mucosa_Subflava_order_list)


# Reorder levels of factor key and y
Mucosa_Subflava_long_df$key <- as.factor(Mucosa_Subflava_long_df$key)
Mucosa_Subflava_long_df$y <- as.factor(Mucosa_Subflava_long_df$y)
Mucosa_Subflava_long_df$key <- factor(Mucosa_Subflava_long_df$key, levels = Mucosa_Subflava_order_list)  
Mucosa_Subflava_long_df$y <- factor(Mucosa_Subflava_long_df$y, levels = Mucosa_Subflava_order_list)  


# rename z to ANI
Mucosa_Subflava_long_df <- Mucosa_Subflava_long_df %>% 
  rename(ANI = "z")

# Then, create the plot:
plot <- ggplot(Mucosa_Subflava_long_df, aes(key, y)) +
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
  theme(axis.text.x = element_text(size = 4, angle = 45, hjust = 1, color = colors_Mucosa_Subflava$colors),
        axis.text.y = element_text(size = 4, color = colors_Mucosa_Subflava$colors),
        axis.ticks=element_blank())

ggsave(file = "/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/ANI/G_0636/G_0636_ANI_heatmap_colors_Mucosa_Subflava_test.pdf", plot = plot, width = 11, height = 8, limitsize = FALSE)

