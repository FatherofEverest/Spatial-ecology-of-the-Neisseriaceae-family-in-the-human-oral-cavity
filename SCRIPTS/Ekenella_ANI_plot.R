# Subset keys dataframe
Eikenella_keys <- keys[538:575, 1, drop = FALSE]

# Convert the dataframe column to a vector
Eikenella_pangenomeIDs <- Eikenella_keys[,1]

# Subset your matrix
Eikenella_long_df <- long_df %>% 
  filter(key %in% Eikenella_pangenomeIDs) %>% 
  filter(y %in% Eikenella_pangenomeIDs) %>% 
  droplevels()


Eikenella_order <- long_df %>% 
  filter(key %in% Eikenella_pangenomeIDs) %>% 
  filter(y %in% Eikenella_pangenomeIDs) %>% 
  droplevels() %>% 
  select(key)

Eikenella_order_list <- levels(Eikenella_order$key)


colors_Eikenella <- as.data.frame(levels(Eikenella_long_df$key))
colors_Eikenella <- colors_Eikenella %>% 
  rename(pangenomeIDs = "levels(Eikenella_long_df$key)") %>% 
  mutate(colors = case_when(grepl('E_glucosivorans', pangenomeIDs) ~ 'magenta',
                            grepl('E_sp', pangenomeIDs) ~ 'darkgreen',
                            grepl('E_exigua', pangenomeIDs) ~ 'red',
                            grepl('E_corrodens', pangenomeIDs) ~ 'midnightblue',
                            grepl('E_halliae', pangenomeIDs) ~ 'chartreuse',
                            grepl('E_sp_str_NML02_A_017', pangenomeIDs) ~ 'orange'))

# remove GCA id from names
colors_Eikenella$pangenomeIDs <- gsub("_id_.*", "", colors_Eikenella$pangenomeIDs)
Eikenella_long_df$key <- gsub("_id_.*", "", Eikenella_long_df$key)
Eikenella_long_df$y <- gsub("_id_.*", "", Eikenella_long_df$y)
Eikenella_long_df$key <- gsub("E_sp_str_NML02_A_017", "E_longinqua_str_NML02_A_017", Eikenella_long_df$key)
Eikenella_long_df$y <- gsub("E_sp_str_NML02_A_017", "E_longinqua_str_NML02_A_017", Eikenella_long_df$y)
Eikenella_order_list <- gsub("_id_.*", "", Eikenella_order_list)
Eikenella_order_list <- gsub("E_sp_str_NML02_A_017", "E_longinqua_str_NML02_A_017", Eikenella_order_list)


# Reorder levels of factor key and y
Eikenella_long_df$key <- as.factor(Eikenella_long_df$key)
Eikenella_long_df$y <- as.factor(Eikenella_long_df$y)
Eikenella_long_df$key <- factor(Eikenella_long_df$key, levels = Eikenella_order_list)  
Eikenella_long_df$y <- factor(Eikenella_long_df$y, levels = Eikenella_order_list)  


# rename z to ANI
Eikenella_long_df <- Eikenella_long_df %>% 
  rename(ANI = "z")

# Then, create the plot:
plot <- ggplot(Eikenella_long_df, aes(key, y)) +
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
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1, color = colors_Eikenella$colors),
        axis.text.y = element_text(size = 10, color = colors_Eikenella$colors),
        axis.ticks=element_blank())

ggsave(file = "/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/ANI/G_0636/G_0636_ANI_heatmap_colors_Eikenella_test.pdf", plot = plot, width = 11, height = 8, limitsize = FALSE)

