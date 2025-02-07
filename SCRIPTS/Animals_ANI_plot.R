
master_ANI_pangenome_ID_order <- as.data.frame(levels(long_df$key))

# Subset keys dataframe
Animal_keys <- master_ANI_pangenome_ID_order[576:6366,, drop = FALSE]


# Convert the dataframe column to a vector
Animal_pangenomeIDs <- Animal_keys[,1]

# Subset your matrix
Animal_long_df <- long_df %>% 
  filter(key %in% Animal_pangenomeIDs) %>% 
  filter(y %in% Animal_pangenomeIDs) %>% 
  droplevels()


Animal_order <- long_df %>% 
  filter(key %in% Animal_pangenomeIDs) %>% 
  filter(y %in% Animal_pangenomeIDs) %>% 
  droplevels() %>% 
  select(key)

Animal_order_list <- levels(Animal_order$key)


colors_Animal <- as.data.frame(levels(Animal_long_df$key))
colors_Animal <- colors_Animal %>% 
  rename(pangenomeIDs = "levels(Animal_long_df$key)") %>% 
  mutate(colors = case_when(grepl('N_sp_str', pangenomeIDs) ~ 'blue',
                            grepl('N_wadsworthii', pangenomeIDs) ~ 'black',
                            grepl('N_yangbaofengii', pangenomeIDs) ~ 'black',
                            grepl('N_musculi', pangenomeIDs) ~ 'black',
                            grepl('N_canis', pangenomeIDs) ~ 'black',
                            grepl('N_zoodegmatis', pangenomeIDs) ~ 'black',
                            grepl('N_dentiae', pangenomeIDs) ~ 'darkseagreen3',
                            grepl('N_dumasiana', pangenomeIDs) ~ 'midnightblue',
                            grepl('N_animaloris', pangenomeIDs) ~ 'black',
                            grepl('N_perflava', pangenomeIDs) ~ 'red',
                            grepl('N_lisongii', pangenomeIDs) ~ 'black',
                            grepl('N_chenwenguii', pangenomeIDs) ~ 'black',
                            grepl('N_animalis', pangenomeIDs) ~ 'black',
                            grepl('N_bacilliformis', pangenomeIDs) ~ 'darkorange1',
                            grepl('N_sp_HMT_020', pangenomeIDs) ~ 'magenta',
                            grepl('K_potus', pangenomeIDs) ~ 'black',
                            grepl('N_shayeganii', pangenomeIDs) ~ 'darkgreen',
                            grepl('K_sp', pangenomeIDs) ~ 'firebrick3',
                            grepl('N_brasiliensis', pangenomeIDs) ~ 'gray60'))


# Deal with duplicate strains N_lisongii_str_ZJ106
#N_lisongii_str_ZJ106_id_GCA_021728675_1 N_lisongii_str_ZJ106_A_id_GCA_021728675_1
#N_lisongii_str_ZJ106_id_GCA_028463985_1 N_lisongii_str_ZJ106_B_id_GCA_028463985_1

colors_Animal$pangenomeIDs <- gsub("N_lisongii_str_ZJ106_id_GCA_021728675_1", "N_lisongii_str_ZJ106_A_id_GCA_021728675_1", colors_Animal$pangenomeIDs)
colors_Animal$pangenomeIDs <- gsub("N_lisongii_str_ZJ106_id_GCA_028463985_1", "N_lisongii_str_ZJ106_B_id_GCA_028463985_1", colors_Animal$pangenomeIDs)
Animal_long_df$key <- gsub("N_lisongii_str_ZJ106_id_GCA_021728675_1", "N_lisongii_str_ZJ106_A_id_GCA_021728675_1", Animal_long_df$key)
Animal_long_df$y <- gsub("N_lisongii_str_ZJ106_id_GCA_021728675_1", "N_lisongii_str_ZJ106_A_id_GCA_021728675_1", Animal_long_df$y)
Animal_long_df$key <- gsub("N_lisongii_str_ZJ106_id_GCA_028463985_1", "N_lisongii_str_ZJ106_B_id_GCA_028463985_1", Animal_long_df$key)
Animal_long_df$y <- gsub("N_lisongii_str_ZJ106_id_GCA_028463985_1", "N_lisongii_str_ZJ106_B_id_GCA_028463985_1", Animal_long_df$y)
Animal_order_list <- gsub("N_lisongii_str_ZJ106_id_GCA_021728675_1", "N_lisongii_str_ZJ106_A_id_GCA_021728675_1", Animal_order_list)
Animal_order_list <- gsub("N_lisongii_str_ZJ106_id_GCA_028463985_1", "N_lisongii_str_ZJ106_B_id_GCA_028463985_1", Animal_order_list)


#N_dentiae_str_DSM_19151
#N_dentiae_str_DSM_19151_id_GCA_002108595_1 N_dentiae_str_DSM_19151_A_id_GCA_002108595_1
#N_dentiae_str_DSM_19151_id_GCA_014055005_1 N_dentiae_str_DSM_19151_B_id_GCA_014055005_1
colors_Animal$pangenomeIDs <- gsub("N_dentiae_str_DSM_19151_id_GCA_002108595_1", "N_dentiae_str_DSM_19151_A_id_GCA_002108595_1", colors_Animal$pangenomeIDs)
colors_Animal$pangenomeIDs <- gsub("N_dentiae_str_DSM_19151_id_GCA_014055005_1", "N_dentiae_str_DSM_19151_B_id_GCA_014055005_1", colors_Animal$pangenomeIDs)
Animal_long_df$key <- gsub("N_dentiae_str_DSM_19151_id_GCA_002108595_1", "N_dentiae_str_DSM_19151_A_id_GCA_002108595_1", Animal_long_df$key)
Animal_long_df$y <- gsub("N_dentiae_str_DSM_19151_id_GCA_002108595_1", "N_dentiae_str_DSM_19151_A_id_GCA_002108595_1", Animal_long_df$y)
Animal_long_df$key <- gsub("N_dentiae_str_DSM_19151_id_GCA_014055005_1", "N_dentiae_str_DSM_19151_B_id_GCA_014055005_1", Animal_long_df$key)
Animal_long_df$y <- gsub("N_dentiae_str_DSM_19151_id_GCA_014055005_1", "N_dentiae_str_DSM_19151_B_id_GCA_014055005_1", Animal_long_df$y)
Animal_order_list <- gsub("N_dentiae_str_DSM_19151_id_GCA_002108595_1", "N_dentiae_str_DSM_19151_A_id_GCA_002108595_1", Animal_order_list)
Animal_order_list <- gsub("N_dentiae_str_DSM_19151_id_GCA_014055005_1", "N_dentiae_str_DSM_19151_B_id_GCA_014055005_1", Animal_order_list)


#N_yangbaofengii_str_ZJ930
#N_yangbaofengii_str_ZJ930_id_GCA_017811745_1 N_yangbaofengii_str_ZJ930_A_id_GCA_017811745_1
#N_yangbaofengii_str_ZJ930_id_GCA_011038775_1 N_yangbaofengii_str_ZJ930_B_id_GCA_011038775_1
colors_Animal$pangenomeIDs <- gsub("N_yangbaofengii_str_ZJ930_id_GCA_017811745_1", "N_yangbaofengii_str_ZJ930_A_id_GCA_017811745_1", colors_Animal$pangenomeIDs)
colors_Animal$pangenomeIDs <- gsub("N_yangbaofengii_str_ZJ930_id_GCA_011038775_1", "N_yangbaofengii_str_ZJ930_B_id_GCA_011038775_1", colors_Animal$pangenomeIDs)
Animal_long_df$key <- gsub("N_yangbaofengii_str_ZJ930_id_GCA_017811745_1", "N_yangbaofengii_str_ZJ930_A_id_GCA_017811745_1", Animal_long_df$key)
Animal_long_df$y <- gsub("N_yangbaofengii_str_ZJ930_id_GCA_017811745_1", "N_yangbaofengii_str_ZJ930_A_id_GCA_017811745_1", Animal_long_df$y)
Animal_long_df$key <- gsub("N_yangbaofengii_str_ZJ930_id_GCA_011038775_1", "N_yangbaofengii_str_ZJ930_B_id_GCA_011038775_1", Animal_long_df$key)
Animal_long_df$y <- gsub("N_yangbaofengii_str_ZJ930_id_GCA_011038775_1", "N_yangbaofengii_str_ZJ930_B_id_GCA_011038775_1", Animal_long_df$y)
Animal_order_list <- gsub("N_yangbaofengii_str_ZJ930_id_GCA_017811745_1", "N_yangbaofengii_str_ZJ930_A_id_GCA_017811745_1", Animal_order_list)
Animal_order_list <- gsub("N_yangbaofengii_str_ZJ930_id_GCA_011038775_1", "N_yangbaofengii_str_ZJ930_B_id_GCA_011038775_1", Animal_order_list)


#N_yangbaofengii_str_ZJ785
#N_yangbaofengii_str_ZJ785_id_GCA_011038745_1 N_yangbaofengii_str_ZJ785_A_id_GCA_011038745_1
#N_yangbaofengii_str_ZJ785_id_GCA_014898075_1 N_yangbaofengii_str_ZJ785_B_id_GCA_014898075_1
colors_Animal$pangenomeIDs <- gsub("N_yangbaofengii_str_ZJ785_id_GCA_011038745_1", "N_yangbaofengii_str_ZJ785_A_id_GCA_011038745_1", colors_Animal$pangenomeIDs)
colors_Animal$pangenomeIDs <- gsub("N_yangbaofengii_str_ZJ785_id_GCA_014898075_1", "N_yangbaofengii_str_ZJ785_B_id_GCA_014898075_1", colors_Animal$pangenomeIDs)

Animal_long_df$key <- gsub("N_yangbaofengii_str_ZJ785_id_GCA_011038745_1", "N_yangbaofengii_str_ZJ785_A_id_GCA_011038745_1", Animal_long_df$key)
Animal_long_df$y <- gsub("N_yangbaofengii_str_ZJ785_id_GCA_011038745_1", "N_yangbaofengii_str_ZJ785_A_id_GCA_011038745_1", Animal_long_df$y)

Animal_long_df$key <- gsub("N_yangbaofengii_str_ZJ785_id_GCA_014898075_1", "N_yangbaofengii_str_ZJ785_B_id_GCA_014898075_1", Animal_long_df$key)
Animal_long_df$y <- gsub("N_yangbaofengii_str_ZJ785_id_GCA_014898075_1", "N_yangbaofengii_str_ZJ785_B_id_GCA_014898075_1", Animal_long_df$y)

Animal_order_list <- gsub("N_yangbaofengii_str_ZJ785_id_GCA_011038745_1", "N_yangbaofengii_str_ZJ785_A_id_GCA_011038745_1", Animal_order_list)
Animal_order_list <- gsub("N_yangbaofengii_str_ZJ785_id_GCA_014898075_1", "N_yangbaofengii_str_ZJ785_B_id_GCA_014898075_1", Animal_order_list)




# remove GCA id from names
colors_Animal$pangenomeIDs <- gsub("_id_.*", "", colors_Animal$pangenomeIDs)
Animal_long_df$key <- gsub("_id_.*", "", Animal_long_df$key)
Animal_long_df$y <- gsub("_id_.*", "", Animal_long_df$y)
Animal_order_list <- gsub("_id_.*", "", Animal_order_list)


# Reorder levels of factor key and y
Animal_long_df$key <- as.factor(Animal_long_df$key)
Animal_long_df$y <- as.factor(Animal_long_df$y)
Animal_long_df$key <- factor(Animal_long_df$key, levels = Animal_order_list)  
Animal_long_df$y <- factor(Animal_long_df$y, levels = Animal_order_list)  


# rename z to ANI
Animal_long_df <- Animal_long_df %>% 
  rename(ANI = "z")

# Then, create the plot:
plot <- ggplot(Animal_long_df, aes(key, y)) +
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
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1, color = colors_Animal$colors),
        axis.text.y = element_text(size = 8, color = colors_Animal$colors),
        axis.ticks=element_blank())

ggsave(file = "/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/ANI/G_0636/G_0636_ANI_heatmap_colors_Animal_test.pdf", plot = plot, width = 11, height = 8, limitsize = FALSE)

