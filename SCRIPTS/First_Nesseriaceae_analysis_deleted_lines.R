
```{r}

library(ggplot2)
library(dplyr)
if(!require(gg.gap)) install.packages("gg.gap")
library(gg.gap)

# load data
df <- read.csv("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/DATA/Neisseriaceae_available_genomes_NCBI_August_2023.csv", header = TRUE)

# summarize counts for genus
summary_df <- df %>% 
  group_by(Genus) %>% 
  summarise(n = n())

# Create a ggplot object
p <- ggplot(summary_df, aes(x = reorder(Genus, -n), y = n)) + 
  geom_bar(stat = "identity", fill = "blue", color = "black") +
  labs(x = "Neisseriaceae genera", y = "Frequency") +
  scale_y_continuous(expand = c(0,0)) + 
  theme_classic() +
  theme(axis.text = element_text(color = "black"))+
  theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1, face = "italic")) 

# clean up plot slightly
p <- p +
  theme(
    panel.background = element_rect(fill = "white"), 
    panel.grid = element_blank(),
    axis.line = element_blank()
  )

# add break in y axis
p <- p %>% 
  gg.gap::gg.gap(
    ylim = c(0, 3800), 
    segments = list(c(200, 3700)),
    tick_width = 50,
    c(0.7,0,0.3)
  )

# save plot
ggsave(file = "/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/FIGURES/Plots/Neisseriacea_genera_available_genomes_NCBI.pdf", plot = p, width = 7, height = 5)

```




##### 2.2 Plot of Distribution of genomes from NCBI


```{r Distribution of genomes from NCBI}

library(ggplot2)
library(dplyr)
library(gg.gap)

# load data
NCBI_unfiltered <- read.csv("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/DATA/Temps/Neisseriaceae_NCBI_metadata_un_filtered.csv", header = TRUE)


NCBI_unfiltered_select_taxa <- NCBI_unfiltered %>% 
  filter(Genus == "Neisseria" | Genus == "Kingella" | Genus == "Eikenella" | Genus == "Simonsiella" | Genus == "NA") %>% 
  filter(Type == "REF") %>% 
  group_by(Genus, Species) %>% 
  summarise(n = n()) 

NCBI_unfiltered_select_taxa$Genus <- as.factor(NCBI_unfiltered_select_taxa$Genus)
NCBI_unfiltered_select_taxa$Species <- as.factor(NCBI_unfiltered_select_taxa$Species)

# set order of genera
NCBI_unfiltered_select_taxa$Genus <- factor(NCBI_unfiltered_select_taxa$Genus, levels = c("Neisseria", "Kingella", "Eikenella", "Simonsiella"))


dodge <- position_dodge(width=0.9)

plot <- ggplot(NCBI_unfiltered_select_taxa, aes(x = reorder(Species, -n), y = n, fill = Genus)) +
  geom_bar(stat = 'identity',position = dodge, color = "black") +
  facet_grid(. ~ Genus, scales = "free", space = "free") +
  scale_y_continuous(expand = c(0, 0)) +
  xlab(NULL)+
  ylab("Frequency") +
  theme_classic() +
  theme(text = element_text(size = 12, color = "black"),
        strip.text.x = element_blank(),
        axis.text.x = element_text(angle = 45,vjust = 1, hjust=1, face = "italic")) 

# clean up plot slightly
plot <- plot +
  theme(
    panel.background = element_rect(fill = "white"), 
    panel.grid = element_blank(),
    axis.line = element_blank()
  )

# add break in y axis
plot <- plot %>% 
  gg.gap::gg.gap(
    ylim = c(0, 2400), 
    segments = list(c(150, 1000)),
    tick_width = 50,
    c(0.7,0,0.3)
  )

ggsave("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/FIGURES/Plots/NCBI_genomes_available.pdf", plot, height = 6, width = 11)

```


Retrieve genbank metadata so that we can link the GCA IDs in the filtered NCBI Neisseriaceae reference genome metadata file with the latest ftp download links for genbank assemblies.
```{bash Retrieve Genbank metadata}
wget -q -P /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/DATA/ https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt

```


Upload filtered Neisseriaceae NCBI genome metadata file; obtain list of GCA assembly IDs; subset assembly_summary_genbank.txt file based on that list to obtain the ftp links.
```{r Get Genbank assembly IDs from NCBI metadata}

library(dplyr)

# load metadata 
metadata_df <- read.csv("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/DATA/Neisseriaceae_NCBI_metadata_filtered.csv", header = TRUE)

Assembly_list <- metadata_df %>% 
  select(Assembly) 

write.table(Assembly_list, "/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/DATA/Neisseriaceae_GCA_assembly_IDs.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

```

```{bash Filter Genbank metadata by selected NCBI genomes}
# set working directory
DIR_DATA=/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/DATA/
  
  # Use grep to filter the rows of file1.txt based on the content of file2.txt using the following command:
  grep -Ff $DIR_DATA/Neisseriaceae_GCA_assembly_IDs.txt $DIR_DATA/assembly_summary_genbank.txt > $DIR_DATA/filtered_assembly_summary_genbank.txt

# get original header
head -n 2 $DIR_DATA/assembly_summary_genbank.txt | tail -n 1 > $DIR_DATA/genbank_header.txt

# concatenate header with filtered genbank metadata
cat $DIR_DATA/genbank_header.txt $DIR_DATA/filtered_assembly_summary_genbank.txt >> $DIR_DATA/final_filtered_assembly_summary_genbank.txt

# convert genbank metadata txt format to csv format
tr '\t' ',' < $DIR_DATA/final_filtered_assembly_summary_genbank.txt > $DIR_DATA/final_filtered_assembly_summary_genbank.csv

```

The NCBI and GenBank filtered Neisseriaceae genome metadata files were then manually merged in excel by first sorting by GenBank assembly ID (A-Z) and then copy and paste. This created a xlsx file that has two tabs - one for the GenBank metadata and one for the NCBI metadata.

14 additional genomes were removed based on Genbank metadat (see supplemental materials: Genome curation).


Now, create a table that summarizes the distribution of genomes across species.
```{r Create summary of genome distribution across species}

library(dplyr)

# load curated metadata 
metadata_df <- read.csv("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/DATA/Neisseriaceae_Final_metadata.csv", header = TRUE)


# Species summary
species_summary <- metadata_df %>% 
  group_by(NCBI_Genus, NCBI_Species) %>% 
  summarise(n = n()) %>% 
  rename(Genus = NCBI_Genus,
         Species = NCBI_Species)

write.table(species_summary, "/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/FIGURES/Tables/species_summary.txt", row.names = FALSE, quote = FALSE, sep = "\t")


metadata_df %>% 
  summarise(n = n())

```


```{bash download genomes}

# upload NCBI metadata file onto the server (Run on local machine)
scp -r /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/DATA/Neisseriaceae_Final_metadata.csv jgiacomini@evol5.mbl.edu:/workspace/jmarkwelchlab/P_0003_Neisseriaceae/DATA/Neisseriaceae_Final_metadata.csv

# copy NCBI metadata for desired genomes (RefSeq)
cat $genomeMetadata | awk 'BEGIN{FS=OFS=","}NR==1{print $0}NR>1{print $0}' > $projectMetadata

# make list with gca path for download; GCA download link is in column 50 in this file
cat $genomeMetadata | awk -F',' -v OFS="\t" 'NR>1{print $50}' | awk 'BEGIN{FS=OFS="/"}{print $0,$NF"_genomic.fna.gz"}' | sed -e 's/"//g' > $downloadNCBI

# Download genomes
clusterize -n 5 -l LOGS/download_genomes.log /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/download_genomes.sh

# Check the number of files in the directory. 
# This should be the number of genomes you expected to be downloaded, which in our case is 488
ls $DIR_NCBI | wc -l 

```

```{bash}

clusterize -n 5 -l LOGS/rename_genomes.log /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/rename_genomes.sh

# check - should be 488
ls $DIR_Assemblies | wc -l

```


```{bash first run}

clusterize -n 5 -l LOGS/reformat_genomes_and_anvio_contigDBs.log /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/reformat_genome_deflines_and_build_contig_databases.sh

# for some reason the last 8 genomes did not complete; so I ran the following script to finish
clusterize -n 5 -l LOGS/build_missing_contig_databases.log /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/build_missing_contig_databases.sh

# check - should be 488
ls $DIR_ContigsDB | wc -l

```

```{bash}

clusterize -n 5 -l LOGS/anvio_hmms_trnas_gtdb.log /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/anvio_hmms_trnas_gtdb.sh

clusterize -n 10 -l LOGS/missing_scg_gtdb_taxonomy.log /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/missing_scg_gtdb_taxonomy.sh

```

```{bash generate contigs stats for each genome using anvio}

mkdir 19_Contig_db_stats

clusterize -n 5 -l LOGS/check_anvio_contigDB_hmms.log /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/check_anvio_contigDB_hmms.sh

cat $mainDIR/$projectID/19_Contig_db_stats/Contig_bacteria_71_summary.txt | wc -l

```

```{bash Rename with HOMD names create ITEMS and contig paths}

clusterize -n 5 -l LOGS/Rename_HOMD_create_ITEMS_Contigs_path_file.log /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/Rename_HOMD_create_ITEMS_Contigs_path_file.sh


# check; should be 488 rows
GENOMES=$mainDIR/$projectID/DATA/$projectID-contig_paths.txt
cat $GENOMES | wc -l 

```


```{r}

df <- read.csv("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/NCBI_DATASETS/02_P_0003_Neisseriaceae.csv", header = TRUE)

write.table(df, "/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/NCBI_DATASETS/02_P_0003_Neisseriaceae.tsv", row.names = FALSE, quote = FALSE)

# Create a function that constructs the desired string
construct_string <- function(row) {
  paste0(substr(row["Genus"], 1, 1), "_", 
         row["Species"], "_str_", 
         row["Strain_abrv"], "_id_", 
         gsub("\\.", "_", row["assembly_accession"]))
}

# Apply this function to each row of the dataframe
result_strings <- apply(df, 1, construct_string)

# Print the result
print(result_strings)



```


##### 3.1 CheckM 1 (Old analysis)

screen -S qrsh_sesh
qrsh -pe allslots 20

module purge
module load miniconda/3
source /bioware/miniconda3/bashrc
conda activate anvio-7.1
module load clusters/barhal
module load jbpc
module load python/3.10.11-20230411

```{bash Genome completeness and contamination using CheckM}

projectID=P_0003_Neisseriaceae
mainDIR=/workspace/jmarkwelchlab
DIR_Data=$mainDIR/$projectID/DATA
nameConversions=$DIR_Data/13_${projectID}-name_conversions.txt
DIR_CheckM=$mainDIR/$projectID/26_CHECKM
DIR_Contigs=$mainDIR/$projectID/03_GENOMES_EDITED

mkdir $DIR_CheckM/BINS_DIR

# copy genomes and change names $rawNCBIid
while IFS= read -r line
do
orginalName=$( echo "$line" | awk -F'\t' '{print $2}')
newName=$( echo "$line" | awk -F'\t' '{print $1}')
cp $DIR_Contigs/$orginalName.fa $DIR_CheckM/BINS_DIR/$newName.fa
done < $nameConversions


binDir=$DIR_CheckM/BINS_DIR
outDIR=$DIR_CheckM/${projectID}-taxonomic
outFile=$outDIR/${projectID}
rank=family
taxon=Neisseriaceae


checkm taxon_set $rank $taxon $outFile.lineage.ms.txt

checkm analyze -x fa -t 100 $outFile.lineage.ms.txt $binDir $outDIR

checkm qa -o 1 -f $outFile.lineage.qa.txt --tab_table -t 100 $outFile.lineage.ms.txt $outDIR

cat $DIR_CheckM/${projectID}-taxonomic/${projectID}.lineage.qa.txt | awk -F "\t" '{print $1,$12, $13}'
```



Check checkM results...

```{bash send checkM results to local}

scp -r jgiacomini@evol5.mbl.edu:/workspace/jmarkwelchlab/P_0003_Neisseriaceae/26_CHECKM/P_0003_Neisseriaceae-taxonomic/P_0003_Neisseriaceae.lineage.qa.txt /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/DATA/CheckM_results_Neisseriaceae.txt

```

********
  Need to remove genomes from the list that have completeness less than 90% and contamination greater than 3
********
  
  We then further filtered the set of reference genomes so that each genome included in the set shared no more than 98% Average Nucleotide Identity (ANI) with any other genome, as well as had a completeness of ≥ 90% and a contamination < 3% estimated by CheckM (56).


Contaminated genomes (>3%)

N_sp_str_ZJ104_id_GCA_021728505_1
N_elongata_str_Nel_M001_id_GCA_018437425_1
K_kingae_str_AA392_id_GCA_030180345_1
K_kingae_str_BB237_id_GCA_030182345_1
K_kingae_str_AA392C5_id_GCA_030180305_1


Incomplete genomes (<90%)

N_weixii_str_10022_id_GCA_002327085_1
E_sp_str_NML03_A_027_id_GCA_001648395_1
E_sp_str_NML01_A_086_id_GCA_001648345_1
*E_exigua_str_PXX_id_GCA_008805035_1*
  E_sp_str_NML070372_id_GCA_001648415_1
E_sp_str_NML96_A_049_id_GCA_001648495_1
E_sp_str_NML97_A_109_id_GCA_001648505_1
E_exigua_str_EI_02_id_GCA_009650935_1
E_sp_str_NML080894_id_GCA_001648425_1
E_sp_str_NML120348_id_GCA_001648435_1
E_sp_str_NML99_0057_id_GCA_001648535_1

I decided to keep the E_exigua_str_PXX_id_GCA_008805035_1 genome to serve as a means of distinguishing from E. corrodens. 


```{bash Remomve incomplete and contaminated genomes from data set}

# copy contigs path file
GENOMES=$mainDIR/$projectID/DATA/$projectID-contig_paths.txt
GENOMES_decontaminated=$mainDIR/$projectID/DATA/$projectID-decontaminated-contig_paths.txt
cp $GENOMES $GENOMES_decontaminated

# remove selected genomes listed above
awk 'FNR==NR{a[$1];next} !($1 in a)' DATA/checkM_decontaminated_incomplete_list.txt $GENOMES > $GENOMES_decontaminated

# make meta data file to add to pangenome
layersADD_decontaminated=DATA/$projectID-decontaminated-add_info.items.txt
ITEMS=${projectID}-add_info.items.txt
head -n1 $ITEMS > $layersADD_decontaminated
for genomeTypeID in `cat $GENOMES_decontaminated`
do
cat $ITEMS | grep "$genomeTypeID" >> $layersADD_decontaminated
done

# write pangenome description txt file for anvi-pan-genome --description flag
PAN_DES_decontaminated=DATA/description_decontaminated_pangenome.txt

echo "This pangenome is of Neisseria, Kingella, Eikenella and Simonsiella species. It's purpose is as an initial pangenome that will provide some information about how genomes cluster, which will be used downstream to filter redunadant genomes or thos ethat cluster with non-human taxa. Genomes with less than 90% completion and greater than 3% contamination estimated by checkM were removed." > $PAN_DES_decontaminated

```


Need to annotate two genomes with bac71 - not sure why this failed on the initial run
```{bash}

#N_animaloris_str_C2015003240_id_GCA_003045095_1
#N_gonorrhoeae_str_NCTC13479_id_GCA_900454155_1

# get genome IDs
#cat DATA/$projectID-decontaminated-contig_paths.txt | grep "N_gonorrhoeae_str_NCTC13479_id_GCA_900454155_1"

DIR_ContigsDB=/workspace/jmarkwelchlab/P_0003_Neisseriaceae/04_CONTIGS_DB
numThreads=5
for genome in G_0211 G_0458
do
contigsDB=$DIR_ContigsDB/${genome}-contigs.db
anvi-run-hmms -c $contigsDB -T $numThreads --just-do-it
done



```

Build decontaminated pangenome and ANI table
```{bash}

mkdir /workspace/jmarkwelchlab/P_0003_Neisseriaceae/09_PANGENOME/P_0003_Neisseriaceae-decontaminated-pangenome

clusterize -n 20 -m jgiacomini@forsyth.org -l LOGS/decontaminated_pangenome.log /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/decontaminated_pangenome.sh

```




