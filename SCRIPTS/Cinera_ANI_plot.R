
master_ANI_pangenome_ID_order <- as.data.frame(levels(long_df$key))

# Subset keys dataframe
Cinerea_keys <- master_ANI_pangenome_ID_order[163:176,, drop = FALSE]


# Convert the dataframe column to a vector
Cinerea_pangenomeIDs <- Cinerea_keys[,1]

# Subset your matrix
Cinerea_long_df <- long_df %>% 
  filter(key %in% Cinerea_pangenomeIDs) %>% 
  filter(y %in% Cinerea_pangenomeIDs) %>% 
  droplevels()


Cinerea_order <- long_df %>% 
  filter(key %in% Cinerea_pangenomeIDs) %>% 
  filter(y %in% Cinerea_pangenomeIDs) %>% 
  droplevels() %>% 
  select(key)

Cinerea_order_list <- levels(Cinerea_order$key)


colors_Cinerea <- as.data.frame(levels(Cinerea_long_df$key))
colors_Cinerea <- colors_Cinerea %>% 
  rename(pangenomeIDs = "levels(Cinerea_long_df$key)") %>% 
  mutate(colors = case_when(grepl('N_sp', pangenomeIDs) ~ 'blue',
                            grepl('N_cinerea', pangenomeIDs) ~ 'red',
                            grepl('N_zalophi', pangenomeIDs) ~ 'magenta'))


# Deal with duplicate strains
# colors_Cinerea$pangenomeIDs <- gsub("N_sp_str_F0314_id_GCA_000090875_1", "N_sp_str_F0314_A_id_GCA_000090875_1", colors_Cinerea$pangenomeIDs)
# colors_Cinerea$pangenomeIDs <- gsub("N_sp_str_F0314_id_GCA_005886145_1", "N_sp_str_F0314_B_id_GCA_005886145_1", colors_Cinerea$pangenomeIDs)
# Cinerea_long_df$key <- gsub("N_sp_str_F0314_id_GCA_000090875_1", "N_sp_str_F0314_A_id_GCA_000090875_1", Cinerea_long_df$key)
# Cinerea_long_df$y <- gsub("N_sp_str_F0314_id_GCA_000090875_1", "N_sp_str_F0314_A_id_GCA_000090875_1", Cinerea_long_df$y)
# Cinerea_long_df$key <- gsub("N_sp_str_F0314_id_GCA_005886145_1", "N_sp_str_F0314_B_id_GCA_005886145_1", Cinerea_long_df$key)
# Cinerea_long_df$y <- gsub("N_sp_str_F0314_id_GCA_005886145_1", "N_sp_str_F0314_B_id_GCA_005886145_1", Cinerea_long_df$y)
# Cinerea_order_list <- gsub("N_sp_str_F0314_id_GCA_000090875_1", "N_sp_str_F0314_A_id_GCA_000090875_1", Cinerea_order_list)
# Cinerea_order_list <- gsub("N_sp_str_F0314_id_GCA_005886145_1", "N_sp_str_F0314_B_id_GCA_005886145_1", Cinerea_order_list)


# remove GCA id from names
colors_Cinerea$pangenomeIDs <- gsub("_id_.*", "", colors_Cinerea$pangenomeIDs)
Cinerea_long_df$key <- gsub("_id_.*", "", Cinerea_long_df$key)
Cinerea_long_df$y <- gsub("_id_.*", "", Cinerea_long_df$y)
Cinerea_order_list <- gsub("_id_.*", "", Cinerea_order_list)


# Reorder levels of factor key and y
Cinerea_long_df$key <- as.factor(Cinerea_long_df$key)
Cinerea_long_df$y <- as.factor(Cinerea_long_df$y)
Cinerea_long_df$key <- factor(Cinerea_long_df$key, levels = Cinerea_order_list)  
Cinerea_long_df$y <- factor(Cinerea_long_df$y, levels = Cinerea_order_list)  


# rename z to ANI
Cinerea_long_df <- Cinerea_long_df %>% 
  rename(ANI = "z")

# Then, create the plot:
plot <- ggplot(Cinerea_long_df, aes(key, y)) +
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
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, color = colors_Cinerea$colors),
        axis.text.y = element_text(size = 12, color = colors_Cinerea$colors),
        axis.ticks=element_blank())

ggsave(file = "/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/ANI/G_0636/G_0636_ANI_heatmap_colors_Cinerea_test.pdf", plot = plot, width = 11, height = 8, limitsize = FALSE)

