
master_ANI_pangenome_ID_order <- as.data.frame(levels(long_df$key))

# Subset keys dataframe
Polysaccharea_keys <- master_ANI_pangenome_ID_order[321:356,, drop = FALSE]


# Convert the dataframe column to a vector
Polysaccharea_pangenomeIDs <- Polysaccharea_keys[,1]

# Subset your matrix
Polysaccharea_long_df <- long_df %>% 
  filter(key %in% Polysaccharea_pangenomeIDs) %>% 
  filter(y %in% Polysaccharea_pangenomeIDs) %>% 
  droplevels()


Polysaccharea_order <- long_df %>% 
  filter(key %in% Polysaccharea_pangenomeIDs) %>% 
  filter(y %in% Polysaccharea_pangenomeIDs) %>% 
  droplevels() %>% 
  select(key)

Polysaccharea_order_list <- levels(Polysaccharea_order$key)


colors_Polysaccharea <- as.data.frame(levels(Polysaccharea_long_df$key))
colors_Polysaccharea <- colors_Polysaccharea %>% 
  rename(pangenomeIDs = "levels(Polysaccharea_long_df$key)") %>% 
  mutate(colors = case_when(grepl('N_sp', pangenomeIDs) ~ 'blue',
                            grepl('N_lactamica', pangenomeIDs) ~ 'red',
                            grepl('N_polysaccharea', pangenomeIDs) ~ 'magenta'))


# Deal with duplicate strains
# colors_Polysaccharea$pangenomeIDs <- gsub("N_sp_str_F0314_id_GCA_000090875_1", "N_sp_str_F0314_A_id_GCA_000090875_1", colors_Polysaccharea$pangenomeIDs)
# colors_Polysaccharea$pangenomeIDs <- gsub("N_sp_str_F0314_id_GCA_005886145_1", "N_sp_str_F0314_B_id_GCA_005886145_1", colors_Polysaccharea$pangenomeIDs)
# Polysaccharea_long_df$key <- gsub("N_sp_str_F0314_id_GCA_000090875_1", "N_sp_str_F0314_A_id_GCA_000090875_1", Polysaccharea_long_df$key)
# Polysaccharea_long_df$y <- gsub("N_sp_str_F0314_id_GCA_000090875_1", "N_sp_str_F0314_A_id_GCA_000090875_1", Polysaccharea_long_df$y)
# Polysaccharea_long_df$key <- gsub("N_sp_str_F0314_id_GCA_005886145_1", "N_sp_str_F0314_B_id_GCA_005886145_1", Polysaccharea_long_df$key)
# Polysaccharea_long_df$y <- gsub("N_sp_str_F0314_id_GCA_005886145_1", "N_sp_str_F0314_B_id_GCA_005886145_1", Polysaccharea_long_df$y)
# Polysaccharea_order_list <- gsub("N_sp_str_F0314_id_GCA_000090875_1", "N_sp_str_F0314_A_id_GCA_000090875_1", Polysaccharea_order_list)
# Polysaccharea_order_list <- gsub("N_sp_str_F0314_id_GCA_005886145_1", "N_sp_str_F0314_B_id_GCA_005886145_1", Polysaccharea_order_list)


# remove GCA id from names
colors_Polysaccharea$pangenomeIDs <- gsub("_id_.*", "", colors_Polysaccharea$pangenomeIDs)
Polysaccharea_long_df$key <- gsub("_id_.*", "", Polysaccharea_long_df$key)
Polysaccharea_long_df$y <- gsub("_id_.*", "", Polysaccharea_long_df$y)
Polysaccharea_order_list <- gsub("_id_.*", "", Polysaccharea_order_list)


# Reorder levels of factor key and y
Polysaccharea_long_df$key <- as.factor(Polysaccharea_long_df$key)
Polysaccharea_long_df$y <- as.factor(Polysaccharea_long_df$y)
Polysaccharea_long_df$key <- factor(Polysaccharea_long_df$key, levels = Polysaccharea_order_list)  
Polysaccharea_long_df$y <- factor(Polysaccharea_long_df$y, levels = Polysaccharea_order_list)  


# rename z to ANI
Polysaccharea_long_df <- Polysaccharea_long_df %>% 
  rename(ANI = "z")

# Then, create the plot:
plot <- ggplot(Polysaccharea_long_df, aes(key, y)) +
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
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, color = colors_Polysaccharea$colors),
        axis.text.y = element_text(size = 12, color = colors_Polysaccharea$colors),
        axis.ticks=element_blank())

ggsave(file = "/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/ANI/G_0636/G_0636_ANI_heatmap_colors_Polysaccharea_test.pdf", plot = plot, width = 11, height = 8, limitsize = FALSE)

