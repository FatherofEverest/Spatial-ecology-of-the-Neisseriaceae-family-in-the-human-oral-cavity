#!/bioware/miniconda3/envs/anvio-7.1/bin/Rscript

## ---------------------------
##
## Script name: Calculating_mean_gene_Shannon_entropy.R
##
## Purpose of script: Calculate mean gene shannon entropy for each genome in each oral site and generate ridgeline plots
## Author: Dr. Jonathan Giacomini
##
## Date Created: 05-12-2023
##
## Copyright (c) Jonathan Giacomini, 2023
## Email: jonjgiacomini@gmail.com
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

## load up the packages we will need: 

library(ggplot2)
# install.packages("ggplot2") installed 5/12/2023 version 3.4.2
library(ggridges) # installed 5/12/2023
#install.packages("ggridges")
library(viridis) # installed 5/12/2023
#install.packages("viridis")
library(dplyr)


## ---------------------------

# The following is an example data frame...
# sample_id = rep(1:2, each = 50)
# corresponding_gene_call = rep(1:10, 10)
# entropy = rnorm(n = 100, mean = 0.5, sd = 0.1)
# df1 = data.frame(sample_id, corresponding_gene_call, entropy)
# colnames(df1) = c("sample_id", "corresponding_gene_call", "entropy")

## ---------------------------

genomes = list("H_parainfluenzae_str_M1C160_1_id_GCA_014931275_1",
             "H_parainfluenzae_str_CCUG_58848_id_GCA_001679405_1",
             "A_aphrophilus_str_C2008003249_id_GCA_003252995_1",
             "A_kilianii_str_PN_528_id_GCA_003130255_1",
             "A_segnis_str_NCTC10977_id_GCA_900476035_1",
             "A_sp_HMT_458_str_W10330_id_GCA_000466335_1",
             "H_haemolyticus_str_60971_B_Hi_3_id_GCA_006439235_1",
             "H_parahaemolyticus_str_C2010039593_id_GCA_003253075_1",
             "H_paraphrohaemolyticus_str_NCTC10671_id_GCA_900451065_1")

sites = list("BM", "TD", "PP")
  
for (genome in genomes) {
  for (site in sites) {
    
    # load data frame $DIR_Variability/${genome}/${site}/${genome}_SNVs.txt
    df1 <- read.table(paste0("/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/27_VARIABILITY/Intra_species_diversity/08_VARIABILITY/",genome,"/",site,"/",genome,"_SNVs.txt"), sep = "\t", header = TRUE)
    
    # subset data frame
    df1 <- df1 %>% 
      select(sample_id, corresponding_gene_call, entropy)

    # remove corresponding_gene_call that is -1, which means position not found in gene
    df1 <- df1 %>% 
     dplyr::filter(corresponding_gene_call != -1)

    # calculate mean shannon entropy for each gene for each sample
    df_final <- df1 %>% 
      dplyr::group_by(sample_id, corresponding_gene_call) %>% 
      dplyr::summarise(MGSE = mean(entropy),
                       MGSE_SD = sd(entropy),
                       N = n())
    df_final$sample_id <- as.factor(df_final$sample_id)
    
    # Reorder y axis by desc peaks so we can overlap them a little (contreol overlap with scale)
    df_final$sample_id <- reorder(df_final$sample_id, desc(df_final$MGSE))
    
    # Plot
    plot <- ggplot(df_final, aes(x = MGSE, y = sample_id, fill = after_stat(x))) +
      geom_density_ridges_gradient(scale = 3, rel_min_height = 0.001, size = 0.4) +
      scale_fill_viridis_c(name = "MGSE", option = "C") +
      scale_x_continuous(expand = c(0, 0), limits = c(-0.2, 1)) +
      scale_y_discrete(expand = expansion(mult = c(0.01, .7))) +
      xlab(NULL) + 
      ylab(NULL) + 
      theme_ridges(center_axis_labels = TRUE) +
      theme(legend.position="none",
            axis.text.y = element_blank()) 
    
    # save plot
    ggsave(paste0("/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/27_VARIABILITY/Intra_species_diversity/08_VARIABILITY/",genome,"/",site,"/ridgeline_plot.pdf"), plot = plot, width = 5, height = 5)
    
    # save test plot
    #ggsave(paste0("/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/27_VARIABILITY/Intra_species_diversity/08_VARIABILITY/test_ridgeline_plot.pdf"), plot = plot, width = 5, height = 10)
    
    # save data frame that summarises MGSE 
    write.table(df_final,paste0("/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/27_VARIABILITY/Intra_species_diversity/08_VARIABILITY/",genome,"/",site,"/MGSE_summary.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
  }
}

