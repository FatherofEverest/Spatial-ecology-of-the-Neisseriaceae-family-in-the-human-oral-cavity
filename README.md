---
title: "Neisseriaceae Metapangenomics in the Human Oral Cavity"
author: "J. J. Giacomini"
date: "2024-04-25"
output: html_document
---

# Load modules on Barhal

The following work was performed on either 1: a computer called Barhal, located at the Marine Biological Laborator in Woods Hole, MA, or 2: a MacBook Pro (14-inch, 2012) Apple M1 Pro with 32 GB RAM.

The specs for Barhal were:

 - 40 64 bit cores / 4 Intel(R) Xeon(R) CPU E5-4610 v3 @ 1.70GHz processors
 - 1763G RAM
 - Rocky Linux release 8.7 (Green Obsidian)
 
 
Within a Barhal screen session, the following modules were run, unless otherwise noted, to create the working environment. 

```{bash, eval=FALSE}

# Loading anvio version 8
module load anvio
module load clusters/barhal
module load jbpc

```

# 1. Genome selection 

To create a reference genome collection representing the diversity of Neisseriaceae in the human oral microbiome, we obtained publicly available RefSeq genomes from the National Center for Biotechnology Information (NCBI) database (downloaded on 2023-09-21) for Neisseria, Eikenella, Kingella and Simonsiella species. We used the NCBI datasets program to retrieve metadata and ftp links used to download the genomes. In total, we downloaded 3937 genomes, of which 2262 were N. meningitidis, 1055 were N. gonorrhoeae, and 620 were other Neisseria, Eikenella, Kingella and Simonsiella taxa. Details about the genomes, including their genus, species, strain, BioSample, BioProject, isolation host, isolation site, RefSeq status, type strain, disease-association, and submitter ID can be found in Supplemental Data: Table S1.

We then performed quality control and dereplication steps to ensure that each genome in the collection had a completeness of at least 90%, had a contamination level below 5% as estimated by CheckM2 (22) (see Supplemental Data: Table S2), and had no more than 98% Average Nucleotide Identity (ANI) with any other genome (see Supplemental Data: Table S3). This process resulted in a set of 213 high-quality reference genomes representing the diversity of Neisseriaceae found in the human microbiome.


##### 1.1 Distribution of Neisseriacea genomes from NCBI

Please see the new NCBI datasets method which is found in the *NCBI_command_line_datasets_tool.Rmd* file

##### 1.2 Set up directories on MBL server

```{bash, eval=FALSE variables for the environment}

# define variables unique to this project 
projectID="P_0003_Neisseriaceae"
mainDIR="/workspace/jmarkwelchlab"
meta=$mainDIR/$projectID/DATA/Neisseriaceae_Final_metadata.csv
ITEMS=${projectID}-add_info.items.txt
seqidHOMD=/workspace/jmarkwelchlab/HOMD_INFO/SEQID_info.txt
DIR_SITE=/storage/data/00_Public_Metagenomes/oral/HMP/
hmpID=/workspace/jmarkwelchlab/HMP_METAGENOMES/METADATA/samples_id-QC.txt
hmpMetadata=/workspace/jmarkwelchlab/HMP_METAGENOMES/METADATA/sequencing-metadata-count-QC.tsv
minContigSIZE=300
TH=15

# directories (variables)
DIR_Data=DATA
DIR_Reads=00_READS
DIR_NCBI=01_NCBI_GENOMES
DIR_MEREN=01_MEREN_GENOMES
DIR_Assemblies=02_ASSEMBLIES
DIR_Contigs=03_GENOMES_EDITED
DIR_ContigsDB=04_CONTIGS_DB
DIR_Mapping=05_MAPPING
DIR_SinglePROF=06_SINGLE_PROFILE
DIR_MergedPROF=07_MERGED_PROFILE
DIR_SummaryPROF=08_PROFILE_SUMMARY
DIR_Pangenome=09_PANGENOME
DIR_SummaryPAN=10_PANGENOME_SUMMARY
DIR_Gene_calls=11_GENE_CALLS
DIR_Annotation=12_FUNCTIONAL_ANNOTATION
DIR_DetectionGENOMES=13_DETECTED_GENOMES
DIR_Phylo=14_PHYLOGENOMICS
DIR_Derep=15_DE_REPLICATION
DIR_var=22_GENE_LEVEL
DIR_CheckM=26_CHECKM

ITEMS
# files (variables)
genomeMetadata=$meta
projectMetadata=$DIR_Data/02_${projectID}.csv
downloadNCBI=$DIR_Data/03_${projectID}-download_url.txt
rawNCBIid=$DIR_Data/04_${projectID}-ncbi_raw_genomes_id.txt
genusList=$DIR_Data/05_${projectID}-genus_list.txt
checkNCBIid=$DIR_Data/06_${projectID}-check_genomes_id.txt
newHOMDid=$DIR_Data/07_${projectID}-homd_id.txt
genomesNCBIid=$DIR_Data/08_${projectID}-ncbi_genomes_id.txt
duplicateStrainNames=$DIR_Data/09_${projectID}-duplicate_strain_names.txt
duplicateStrainInfo=$DIR_Data/10_${projectID}-duplicate_strain_names_info.txt
removeDuplicate=$DIR_Data/11_${projectID}-duplicate_strain_removed.txt
magsID=$DIR_Data/12_mags_id-meren.txt
nameConversions=$DIR_Data/13_${projectID}-name_conversions.txt
genomesID=$DIR_Data/id_genomes.txt
genomesALL=$DIR_Data/id_genomes-ALL.txt
genomesRefSeq=$DIR_Data/id_genomes-RefSeq.txt
RAW_ASSEMBLY=$DIR_Assemblies/${projectID}-RAW.fa
PROJECT_CONTIGS=$DIR_Contigs/${projectID}.fa
PROJECT_REPORT=$DIR_Contigs/${projectID}.report.tsv
BINNING=$DIR_Contigs/${projectID}.binning.tsv
DECOMPOSE=$DIR_Contigs/${projectID}.decompose.tsv
CONTIGS_DB=$DIR_ContigsDB/${projectID}-contigs.db
GENE_CALLS=$DIR_Gene_calls/${projectID}-gene_calls.fa
PAN_LAYERS=$DIR_Pangenome/${projectID}-add_info.layers.tsv
samplesMetadata=DATA/${projectID}-samples_metadata.txt
LAYER_ORDERS=DATA/${projectID}-layer_orders.txt

genomesALL=DATA/id_genomes-ALL.txt
nameConversions=$DIR_Data/13_${projectID}-name_conversions.txt
iDir=$DIR_Pangenome/internal_annotated_${projectID}
genomeDB=$iDir/${projectID}-GENOMES.db
panDB_detection=$iDir/${projectID}-RESULTS/${projectID}-PAN_detection.db
genomes_98ANI=DATA/id_genomes-98ANI.txt

PanProject=${projectID}-98ANI
PanDir=$DIR_Pangenome/$PanProject
GENOMES_98ANI_DEREP=$PanDir/${PanProject}.txt
GENOMES_98ANI_DEREP_DB=$PanDir/${PanProject}-GENOMES.db
PAN_DIR=$PanDir/${PanProject}-RESULTS
PAN_DB=$PanDir/${PanProject}-RESULTS/${PanProject}-PAN.db
ANI_DIR=$PanDir/${PanProject}-RESULTS/ANI_RESULTS
layersADD_98ANI=$PanDir/${PanProject}-layers_98ANI.tsv
PAN_DES=DATA/description_98ANI_derep_pangenome.txt

```

```{bash, eval=FALSE}
#
mkdir $mainDIR/$projectID && cd $mainDIR/$projectID

# create directories
mkdir $DIR_Data $DIR_Reads $DIR_NCBI $DIR_MEREN $DIR_Assemblies $DIR_Contigs $DIR_ContigsDB $DIR_Mapping $DIR_SinglePROF $DIR_MergedPROF $DIR_SummaryPROF $DIR_Pangenome $DIR_SummaryPAN $DIR_DetectionGENOMES $DIR_Gene_calls $DIR_Annotation $DIR_Phylo $DIR_Derep $DIR_CheckM

# copy samples id and metadata-count
#cp $hmpID samples_id-QC.txt
#cp $hmpMetadata sequencing-metadata-count-QC.tsv
```


##### 1.3 Download genomes

I created a metadata file with information about the 620 genomes that we need (excluding Men and Gon genomes that we are processing separately. The file is called RefSeq_non_men_or_gon.csv which is the second tab of the 10_13_2023_FINAL_merged_datasets_and_genbank_metadata.xlsx file. The RefSeq_non_men_or_gon.csv can be used to download the new genomes.

```{bash, eval=FALSE download new genomes}
# upload NCBI metadata file onto the server (Run on local machine)
scp -r /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/NCBI_DATASETS/RefSeq_non_men_or_gon.csv BARHAL_SERVER:/workspace/jmarkwelchlab/P_0003_Neisseriaceae/DATA/RefSeq_non_men_or_gon_new_genomes.csv

scp -r /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/NCBI_DATASETS/ftp_paths_parents.txt BARHAL_SERVER:/workspace/jmarkwelchlab/P_0003_Neisseriaceae/DATA/new_genomes_ftp_paths_parents.txt

genomeMetadata=/workspace/jmarkwelchlab/P_0003_Neisseriaceae/DATA/RefSeq_non_men_or_gon_new_genomes.csv

# copy NCBI metadata for desired genomes (RefSeq)
cat $genomeMetadata | awk 'BEGIN{FS=OFS=","}NR==1{print $0}NR>1{print $0}' > $projectMetadata

# make url file
cat DATA/new_genomes_ftp_paths_parents.txt | awk 'BEGIN{FS=OFS="/"}{print $0,$NF"_genomic.fna.gz"}' > $downloadNCBI

# Download genomes
clusterize -n 5 -l LOGS/New_genomes_10_13_2023_download_genomes.log /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/download_genomes.sh

```

##### 1.4 Rename genomes

Anvio doesn't like complicated genome IDs so we have to create a list of unique IDs containing a number and the assembly ID separated by "_id_"; use GCA assembly ID column 9

```{bash, eval=FALSE rename new genomes}

projectID="P_0003_Neisseriaceae"
mainDIR="/workspace/jmarkwelchlab"
DIR_NCBI=$mainDIR/$projectID/01_NCBI_GENOMES
DIR_Assemblies=$mainDIR/$projectID/02_ASSEMBLIES
DIR_Data=$mainDIR/$projectID/DATA
projectMetadata=$DIR_Data/02_${projectID}.csv
genomesNCBIid=$DIR_Data/new_08_P_0003_Neisseriaceae-ncbi_genomes_id.txt
rawNCBIid=$DIR_Data/new_04_P_0003_Neisseriaceae-ncbi_raw_genomes_id.txt
genomesID=$DIR_Data/new_id_genomes.txt

# create a list of unique IDs containing a number and the assembly ID separated by "_id_"; use GCA assembly ID column 7
# (e.g. G_0001_id_GCA_000000000.0)
cat $projectMetadata | awk -F',' -v OFS="\t" 'NR>1{ printf "G_""%04i_id_%s\n", NR-1,$1 }' | sed -e 's/"//g'  > $rawNCBIid

# copy
cp $rawNCBIid $genomesNCBIid
cp $rawNCBIid $genomesID

# create new fasta files for each genome and 
# change names of original fasta file to number a number using the $genomesID file
# (e.g. GCA_901873365.1_Aggregatibacter_sp._BgEED05_genomic.fna -> G_0001-RAW.fa)

while IFS= read -r genomes_id
do
gca_id=$(echo $genomes_id | awk -F'_id_' '{print $2}')
old_name=$(find $DIR_NCBI -name "$gca_id*")
NEW_NAME=$(echo "$genomes_id" | awk -F'_id_' -v new_dir="$DIR_Assemblies" '{print new_dir"/"$1"-RAW.fa"}' )
cp $old_name $NEW_NAME
done < $genomesID

# remove GCA number from $genomesID file to match new fasta file names
sed -i 's/_id_.*//' $genomesID

```

##### 1.5 Reformat genome deflines & build anvio contigs databases for each genome

Using ANVIO (anvi-script-reformat-fasta) to change any non-canonical letter with N this is done in batches of 20 in a background process

Anvi'o used Prodigal (v2.6.3) by Hyatt et al (doi:10.1186/1471-2105-11-119) to identify open reading frames in the data. 

```{bash, eval=FALSE new genomes}

clusterize -n 5 -l LOGS/New_620_genomes_reformat_genomes_and_anvio_contigDBs.log /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/reformat_genome_deflines_and_build_contig_databases.sh

```


Looks like one contigs DB failed... 
```{bash check contig databases} 

# check: the number should be 620, but is actually 619
cat LOGS/New_620_genomes_reformat_genomes_and_anvio_contigDBs.log | grep "✓ anvi-gen-contigs-database" | wc -l

# Make a file that shows if each genome has a successful contig db built
cat LOGS/New_620_genomes_reformat_genomes_and_anvio_contigDBs.log | grep "/workspace/jmarkwelchlab/P_0003_Neisseriaceae/04_CONTIGS_DB/G_" | sed -e 's!/workspace/jmarkwelchlab/P_0003_Neisseriaceae/04_CONTIGS_DB/G_!!g' | sed -e 's/-contigs.db,//g' > contigs_check.txt

# remove empty spaces
sed -i 's/[[:space:]]//g' contigs_check.txt

# The following code will tell you which row is missing in the contigs_check.txt
# Create an array with the expected sequence
expected_sequence=($(seq -w 0001 0620))

# Read the lines from your file into an array
readarray -t file_lines < contigs_check.txt

# Iterate through the expected sequence
for i in "${expected_sequence[@]}"; do
    # Check if the current number is in the file lines
    if [[ ! " ${file_lines[*]} " =~ " $i " ]]; then
        echo "Missing row: $i"
    fi
done


#Missing row: 0243


# Now find the error info for that contig. Looks like a disk error...lame!
grep -A 50 "Name .........................................: G_0243" LOGS/New_620_genomes_reformat_genomes_and_anvio_contigDBs.log 

# check db info since there is a file for it
anvi-db-info 04_CONTIGS_DB/G_0243-contigs.db
# Says configure error!

# rebuild the contig db
anvi-gen-contigs-database -f 03_GENOMES_EDITED/${genome}.fa -n G_0243 -o 04_CONTIGS_DB/G_0243-contigs.db -T 5 

# Had to migrate to version 21 for anvio 8
clusterize -n 5 -l LOGS/New_620_genomes_migrate_v21.log /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/migrate_contigs_db_v21.sh 

```



##### 1.6 Run hmms to annotate with marker gene DBs (Bac_71, rRNAs Arch and Bac); tRNAs; GTDB taxonomy

```{bash, eval=FALSE new genomes}

clusterize -n 5 -l LOGS/new_genomes_anvio_hmms_trnas_gtdb.log /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/anvio_hmms_trnas_gtdb.sh

cat LOGS/new_genomes_anvio_hmms_trnas_gtdb.log | grep -A 5 "Done with Bacteria_71 🎊" | grep "Number of hits in annotation dict" | wc -l
# Shows 620 completed, as expected
```


##### 1.7 Contigs stats for each genome using anvio.

This is an important step to make sure that each and every contig has the same hmms annotations. Otherwise there will be issues downstream. 

```{bash, eval=FALSE new genomes}

mkdir 19_Contig_db_stats/OLD
mv 19_Contig_db_stats/*stats 19_Contig_db_stats/OLD/
mv 19_Contig_db_stats/Contig_bacteria_71_summary.txt 19_Contig_db_stats/OLD/

clusterize -n 1 -l LOGS/new_genomes_check_anvio_contigDB_hmms.log /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/check_anvio_contigDB_hmms.sh

```

##### 1.8 HOMD names, ITEMS, Contig paths file

Import HOMD names for genomes, create and "ITEMS" file that will be used to add metadata to the pangenome and create a Contigs path file that is needed for profiling mapped reads downstream.

```{bash new genomes }

# Upload metadata file from local compupter to Barhal
scp -r /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/NCBI_DATASETS/02_P_0003_Neisseriaceae.csv BARHAL_SERVER:/workspace/jmarkwelchlab/P_0003_Neisseriaceae/DATA/02_P_0003_Neisseriaceae.csv

# run script
clusterize -n 5 -l LOGS/new_genomes_Rename_HOMD_create_ITEMS_Contigs_path_file.log /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/Rename_HOMD_create_ITEMS_Contigs_path_file.sh

```

##### 1.9 Send pangenome IDs to local for adding to metadataa

```{bash, eval=FALSE}

scp -r BARHAL_SERVER:/workspace/jmarkwelchlab/P_0003_Neisseriaceae/DATA/13_P_0003_Neisseriaceae-name_conversions.txt /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/DATA/13_P_0003_Neisseriaceae-name_conversions.txt

```


##### 1.10 CheckM 2 

```{bash, eval=FALSE}

# Barhal set-up
screen -S qrsh_sesh_2
qrsh -pe allslots 20

module load checkm2
module load clusters/barhal
module load jbpc


projectID="P_0003_Neisseriaceae"
mainDIR="/workspace/jmarkwelchlab"
DIR_Data=$mainDIR/$projectID/DATA
nameConversions=$DIR_Data/13_${projectID}-name_conversions.txt
DIR_CheckM=$mainDIR/$projectID/26_CHECKM
DIR_Contigs=$mainDIR/$projectID/03_GENOMES_EDITED

mkdir $DIR_CheckM/BINS_DIR

# copy genomes and change names
while IFS= read -r line
do
orginalName=$( echo "$line" | awk -F'\t' '{print $2}')
newName=$( echo "$line" | awk -F'\t' '{print $1}')
cp $DIR_Contigs/$orginalName.fa $DIR_CheckM/BINS_DIR/$newName.fa
done < $nameConversions


checkm2 predict --threads 5  -x .fa --input $mainDIR/$projectID/26_CHECKM/BINS_DIR/ --output-directory $mainDIR/$projectID/26_CHECKM/CHECKM2/Genomes_620


```

Send results to local machine 
```{bash, eval=FALSE Send CheckM 2 results to local machine }

scp -r BARHAL_SERVER:/workspace/jmarkwelchlab/P_0003_Neisseriaceae/26_CHECKM/CHECKM2/Genomes_620 /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/26_CHECKM/CHECKM2/
```


The following genomes were removed from the data set due to low completeness (< 90%) and high contamination (> 5%) as estimated by CheckM 2:

Genome_ID   Completeness    Contamination
K_potus_str_3_SID_1128_id_GCA_022870985_1   80.73   4.29
N_weixii_str_10022_id_GCA_002327085_1   85.22   2.8
K_kingae_str_AA392_id_GCA_030180345_1   89.92   4.41
N_lisongii_str_ZJ104_id_GCA_021728505_1   100   9.4
N_elongata_str_Nel_M001_id_GCA_018437425_1    99.61   8.77

The follwoing code used to remove incomplete and contaminated genomes from 620 genome data set.

```{bash, eval=FALSE}

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

G_0620_name_conv=DATA/13_P_0003_Neisseriaceae-name_conversions.txt
G_0615_name_conv=DATA/G_0615-name_conversions.txt
# remove selected genomes listed above
awk 'FNR==NR{a[$1];next} !($1 in a)' DATA/checkM_decontaminated_incomplete_list.txt $G_0620_name_conv > $G_0615_name_conv


# write pangenome description txt file for anvi-pan-genome --description flag
PAN_DES_decontaminated=DATA/description_decontaminated_pangenome.txt

echo "This pangenome is of Neisseria, Kingella, Eikenella and Simonsiella species. It's purpose is as an initial pangenome that will provide some information about how genomes cluster, which will be used downstream to filter redunadant genomes or thos ethat cluster with non-human taxa. Genomes with less than 90% completion and greater than 3% contamination estimated by checkM were removed." > $PAN_DES_decontaminated

```


##### 1.11. Meningitidis & Gonorrhoeae

Neisseeria meningitidis & gonorrhoeae gemnomes were processed seperately from the others simply due to the extremely large number of publicaly available genomes for these two species. 

###### 1.11.1 Dowload and process genomes

*See file meningitidis_gonorrhoeae_NCBI_genome_dereplication.Rmd*

###### 1.11.2 CheckM2 on Men/Gno genomes

module load checkm2
module load clusters/barhal
module load jbpc

```{bash, eval=FALSE}
projectID="P_0003_Neisseriaceae"
mainDIR="/workspace/jmarkwelchlab"

clusterize -n 20 -m jgiacomini@forsyth.org -l LOGS/checkm2_initial_run.log checkm2 predict --threads 20  -x .fa --input $mainDIR/$projectID/26_CHECKM/BINS_DIR/ --output-directory $mainDIR/$projectID/26_CHECKM/CHECKM2/

mkdir $mainDIR/$projectID/26_CHECKM/CHECKM2/MEN_GNO
clusterize -n 20 -m jgiacomini@forsyth.org -l LOGS/checkm2_men_gno.log checkm2 predict --threads 20  -x .fa --input $mainDIR/$projectID/03_GENOMES_EDITED/MEN_GNO_GENOMES --output-directory $mainDIR/$projectID/26_CHECKM/CHECKM2/MEN_GNO

```


```{bash send checkm2 reuslts to local}

mkdir  /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/26_CHECKM  /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/26_CHECKM/CHECKM2/

scp -r BARHAL_SERVER:/workspace/jmarkwelchlab/P_0003_Neisseriaceae/26_CHECKM/CHECKM2/quality_report.tsv /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/26_CHECKM/CHECKM2/quality_report.tsv


scp -r BARHAL_SERVER:/workspace/jmarkwelchlab/P_0003_Neisseriaceae/26_CHECKM/CHECKM2/MEN_GNO/quality_report.tsv /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/26_CHECKM/CHECKM2/MEN_GNO_quality_report.tsv


scp -r BARHAL_SERVER:/workspace/jmarkwelchlab/P_0003_Neisseriaceae/DATA/meningitidis_gonorrhoeae-name_conversions.txt /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/DATA/meningitidis_gonorrhoeae-name_conversions.txt
```

###### 1.11.3 Dereplicate Men/Gno genomes 

De-replicate of Neisseria meningitidis and Neisseria gonorrhoeae genomes.

Please see file "meningitidis_gonorrhoeae_NCBI_genome_dereplication.Rmd" for explanation of the code used to download and process the genomes. 

At this point all 3317 genomes were successfully downloaded, unzipped, reformatted and contigs databases were generated for each. The next step is dereplication using Anvio and FastANI.

```{bash, eval=FALSE}
anvi-dereplicate-genomes -e $MEN_GNO_GENOMES -o $DIR_Derep/meningitidis_gonorrhoeae_ANI_98_dereplication_fastANI --skip-fasta-report --program fastANI --similarity-threshold 0.98 --cluster-method simple_greedy --representative-method centrality --num-threads 20 --log-file $mainDIR/$projectID/LOGS/meningitidis_gonorrhoeae-anvi-dereplicate-fastANI_98.log

#Check
cat $mainDIR/$projectID/LOGS/meningitidis_gonorrhoeae-anvi-dereplicate-fastANI_98.log | grep "query" | wc -l
```

Run mode .....................................: fastANI

CITATION

Anvi'o will use 'fastANI' by Jain et al. (DOI: 10.1038/s41467-018-07641-9) to
compute ANI. If you publish your findings, please do not forget to properly
credit their work.

[fastANI] Kmer size ..........................: 16
[fastANI] Fragment length ....................: 3,000
[fastANI] Min fraction of alignment ..........: 0.25
[fastANI] Num threads to use .................: 20
[fastANI] Log file path ......................: /workspace/jmarkwelchlab/P_0003_Neisseriaceae/LOGS/meningitidis_gonorrhoeae-anvi-dereplicate-fastANI_98.log

[10 Sep 23 23:21:49 fastANI] Many to Many ...

fastANI similarity metric ....................: calculated
Number of genomes considered .................: 3,317
Number of redundant genomes ..................: 3,296
Final number of dereplicated genomes .........: 21

ANI RESULTS

* Matrix and clustering of 'ani' written to output directory
* Matrix and clustering of 'alignment fraction' written to output directory
* Matrix and clustering of 'mapping fragments' written to output directory
* Matrix and clustering of 'total fragments' written to output directory


There were 3,317 genomes (1055 for N. gonorrhoeae and 2262 for N. meningitidis). I used FastANI, with default settings, via the anvio program anvi-dereplicate-genomes to dereplicate based on a 98% ANI similarity threshold. Using pyANI instead of fastANI would have taken significantly longer. The results are a major reduction in the number of genomes - 3,296 were considered redundant, resulting in a total of 21 clusters. There are 20 clusters for N. meningitidis, one of which contains 2094 genomes, and there is only a single cluster for all 1055 N. gonorrhoeae genomes.  
 
Overall, I am happy to see such a major reduction in the number of genomes, but I wonder if it is any concern that there is only one N. gonorrhoeae cluster. This is not an oral species, so it likely does not matter for mapping. And, more importantly, the Neisseria, Eikenella, Kingella and Simonsiella phylogenies that I built previously using 26 complete N. gonorrhoeae genomes show a monophyletic group for both N. gonorrhoeae. Additionally, the metaphlan results show zero abundance of N. gonorrhoeae for all samples.

As to which genome to pick -- in Anthony's work with Strep, we biased the selection toward (1) closed genomes (complete genomes in NCBI), (2) the type strain, and (3) genomes that can be bought from culture collections, even if not the type strain.

There is one genome in RefSeq assembled from type strain material (N_gonorrhoeae_str_DSM_9188_id_GCA_003315235_1). It is not considered to be a complete genome in NCBI because it has 74 scaffolds with a few gaps throughout. There are a small number of complete genomes that come from the NCTC culture collection (N = 4). Of those, CheckM2 estimates 100% completion for one (N_gonorrhoeae_str_NCTC13799_id_GCA_900186935_1 ) and 99.9% for the other three. All four were estimated to have 0.11% contamination. 
For meningitidis, there are two RefSeq genomes assembled from type strain material, both of which are considered complete by NCBI (N_meningitidis_str_NCTC10025_id_GCA_900638555_1 and N_meningitidis_str_PartJ_Nmeningitidis_RM8376_id_GCA_022869645_1). 


```{bash, eval=FALSE}

# print number of clusters, number of genomes per cluster, and chosen genomes
cat 15_DE_REPLICATION/meningitidis_gonorrhoeae_ANI_98_dereplication_fastANI/CLUSTER_REPORT.txt | awk '{print $1, $2, $3}'

# Make sure there are no other N. gonorrhoeae genomes nested with N meningitidis clusters
cat 15_DE_REPLICATION/meningitidis_gonorrhoeae_ANI_98_dereplication_fastANI/CLUSTER_REPORT.txt | awk 'NR>3' | grep "gonorrhoeae" | wc -l

# Where is the type strain genome for N_meningitidis? Its in cluster_000001
cat 15_DE_REPLICATION/meningitidis_gonorrhoeae_ANI_98_dereplication_fastANI/CLUSTER_REPORT.txt | grep "N_meningitidis_str_NCTC10025_id_GCA_900638555_1" | awk '{print $1}'

# The other type strain genome is also in cluster_000001
cat 15_DE_REPLICATION/meningitidis_gonorrhoeae_ANI_98_dereplication_fastANI/CLUSTER_REPORT.txt | grep "N_meningitidis_str_PartJ_Nmeningitidis_RM8376_id_GCA_022869645_1" | awk '{print $1}'


# send results to local machine for inclusion into metadata file. 
scp -r BARHAL_SERVER:/workspace/jmarkwelchlab/P_0003_Neisseriaceae/15_DE_REPLICATION/meningitidis_gonorrhoeae_ANI_98_dereplication_fastANI/CLUSTER_REPORT.txt /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/15_DE_REPLICATION/meningitidis_gonorrhoeae_ANI_98_dereplication_fastANI-CLUSTER_REPORT.txt

scp -r BARHAL_SERVER:/workspace/jmarkwelchlab/P_0003_Neisseriaceae/15_DE_REPLICATION/meningitidis_gonorrhoeae_ANI_98_dereplication_fastANI/SIMILARITY_SCORES/fastANI_ani.txt /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/15_DE_REPLICATION/meningitidis_gonorrhoeae_ANI_98_dereplication_fastANI_ani.txt


mkdir /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/ANI/MEN_GON_SIMILARITY_SCORES
scp -r BARHAL_SERVER:/workspace/jmarkwelchlab/P_0003_Neisseriaceae/15_DE_REPLICATION/meningitidis_gonorrhoeae_ANI_98_dereplication_fastANI/SIMILARITY_SCORES/ /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/ANI/MEN_GON_SIMILARITY_SCORES/

```

The above information satisfies how to handle the two large clusters, but we need a strategy to sort through the other N meningitidis clusters to pick genomes by (1) closed genomes (complete genomes in NCBI), (2) the type strain, and (3) genomes that can be bought from culture collections.For the strain that we would select, we should also inspect how it fares by the average similarity metric.

I first generated a txt file for each cluster that contain a comma seperated list of genomes included in the cluster. I then transformed the comma-separated elements in the cluster_00000X.txt files into a column; We can use the tr command to replace commas with newline characters.
```{bash, eval=FALSE}

cd /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/15_DE_REPLICATION/Men_Gon_clusters/

for cluster in 000001 000002 000003 000004 000005 000006 000008 000009 000016 000018
do

# First extract the list of genomes per each cluster
cat /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/15_DE_REPLICATION/meningitidis_gonorrhoeae_ANI_98_dereplication_fastANI-CLUSTER_REPORT.txt | grep "cluster_$cluster"  |awk -F"\t" '{print$4}' > cluster_$cluster.txt

# Then transform the comma-separated elements in the cluster_00000X.txt files into a column;
# We can use the tr command to replace commas with newline characters.
cat cluster_$cluster.txt | tr ',' '\n' > transformed_cluster_$cluster.txt
done


```


N. meningitidis and N. gonorrhoeae ANI heatmap...

```{r}
# packages
library(ggplot2)
library(tidyr)
library(reshape2)
library(phytools)
library(tibble)
library(dendextend)
library(ape)
library("dplyr")

# load ANI data
df <-read.table("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/ANI/MEN_GON_SIMILARITY_SCORES/fastANI_ani.txt", header = TRUE)

# load ANI tree
df_newick <-  ape::read.tree("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/ANI/MEN_GON_SIMILARITY_SCORES/fastANI_ani.newick")

# order ANI data by ANI tree
df_newick_orders <- phytools::compute.mr(df_newick, type = "matrix")
df_newick_rownames <- rev(rownames(df_newick_orders))
data_ordered <- df[ order(match(df$key, df_newick_rownames)), ]
df_newick_rownames_df <- as.data.frame(df_newick_rownames)
cols_A<- ncol(df) -1
df_newick_rownames_df2 <- as.data.frame(matrix(df_newick_rownames_df$df_newick_rownames, ncol = cols_A, byrow = TRUE))
names(df_newick_rownames_df2) <- df_newick_rownames_df2[1,]
df_newick_rownames_df3 <- df_newick_rownames_df2[-1,]
df_newick_rownames_df3 <- df_newick_rownames_df3 %>% 
  add_column(key = NA, .before = 1)
data_ordered2<-data_ordered[names(df_newick_rownames_df3)]

cols2<- ncol(data_ordered2) 
long_df <- reshape2::melt(data_ordered2,
                          id.vars=c("key"),
                          measure.vars=colnames(data_ordered2[2:cols2]),
                          variable.name="y",
                          value.name="z")
mylevels1 <- df_newick$tip.label
long_df$key <- factor(long_df$key,levels=mylevels1)
long_df$y <- factor(long_df$y, levels=mylevels1)

long_df$rounded_z <- round(long_df$z, digits = 2)

# plot
plot<-ggplot(long_df, aes(key,y)) +
  geom_tile(aes(fill = rounded_z)) + 
  scale_fill_gradient2(low = "darkblue",
                       mid = "white",
                       high = "darkred",
                       midpoint = 0.98,
                       limits =c(0.97, 1),
                       na.value="darkblue")+
  #geom_text(aes(label = format(round(z, digits=3), nsmall = 3)),size=2.75) +
  ylab("") +  
  xlab("") + 
  theme_classic() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x= element_blank()) 


ggsave(file = "/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/ANI/MEN_GON_SIMILARITY_SCORES/Small_MEN_GON_ANI_heatmap_rounded.pdf", plot = plot, width = 11, height = 8)
ggsave(file = "/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/ANI/MEN_GON_SIMILARITY_SCORES/Small_MEN_GON_ANI_heatmap.pdf", plot = plot, width = 11, height = 8)


# rows <- as.data.frame(df_newick$tip.label)
# height=nrow(rows) *0.15
# width=height*1.2
#ggsave(file = "/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/ANI/MEN_GON_SIMILARITY_SCORES/MEN_GON_ANI_heatmap.pdf", plot = plot, width = width, height = height, limitsize = FALSE)


```

Now merge each list with the metadata file, keeping only the rows that match the GCA IDs (need to extract) in the list
```{r}
# packages
library(ggplot2)
library(tidyr)
library(reshape2)
library(phytools)
library(tibble)
library(dendextend)
library(ape)
library("dplyr")

metadata <- read.csv("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/DATA/10_13_2023_FINAL_merged_datasets_and_genbank_metadata_tab.csv", header = TRUE)
  
ani <- read.table("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/15_DE_REPLICATION/meningitidis_gonorrhoeae_ANI_98_dereplication_fastANI_ani.txt", header = TRUE, sep = "\t")

clusterIDs <- list("000001", "000002", "000003", "000004", "000005", "000006", "000008", "000009", "000016", "000018")

for (clusterID in clusterIDs) {
  
    list <- read.table(paste0("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/15_DE_REPLICATION/Men_Gon_clusters/transformed_cluster_",clusterID,".txt"), header = FALSE)
    
    # For each pangenome ID in the list, extract the GCA ID
    list <- list %>% 
      mutate(assembly_accession = gsub(".*_id_", "", V1))
    
    list$assembly_accession <- gsub("_1", ".1", list$assembly_accession)
    list$assembly_accession <- gsub("_2", ".2", list$assembly_accession)
    
    # Now merge with the metadata file based on assembly_accession
    merged_list <- merge(metadata, list, by = "assembly_accession")
    
    # save
    write.csv(merged_list, paste0("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/15_DE_REPLICATION/Men_Gon_clusters/transformed_cluster_",clusterID,"_with_metadata.csv"), row.names = FALSE, quote = FALSE)
    
    # Also, for each cluster, subset the master ANI matrix
    ani_filtered <- ani %>% 
      filter(key %in% list$V1)
      
    # Extract the list of pangenome IDs from the 'list' dataframe
    selected_columns <- as.character(list$V1)
    
    # Make sure to include the "key" column as well, since you want to retain it
    selected_columns <- c("key", selected_columns)
    
    # Subset the 'ani' dataframe using the selected columns
    ani_filtered_2 <- ani_filtered[, selected_columns, drop = FALSE]
    
    # save ani results
    write.table(ani_filtered_2, paste0("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/15_DE_REPLICATION/Men_Gon_clusters/transformed_cluster_",clusterID,"_ANI_matrix.txt"), row.names = FALSE, quote = FALSE, sep = "\t")

    
    wide_df2 <- ani_filtered_2[,-1]
    rownames(wide_df2) <- ani_filtered_2[,1]

    # convert to matrix for heirachrical clustering
    wide_matrix <- as.matrix(wide_df2)
    
    # calculate Bray-Curtis distance among samples
    dist_df <- vegan::vegdist(wide_matrix, method = "bray")
    
    #Replace na values with 0 using is.na()
    dist_df[dist_df=="NaN"]<-0
    
    # Hierarchical clustering using Complete Linkage
    hc <- hclust(dist_df, method = "complete" )
      
    # transform to longform 
    cols2<- ncol(ani_filtered_2) 
    long_df <- reshape2::melt(ani_filtered_2,
                              id.vars=c("key"),
                              measure.vars=colnames(ani_filtered_2[2:cols2]),
                              variable.name="y",
                              value.name="z")
    
    
    # reorder genomes for triples plot based on hclust
    genomes_order <- hc$labels[hc$order]
    long_df$key <- factor(long_df$key, levels = genomes_order)  
    long_df$y <- factor(long_df$y, levels = genomes_order)  
    
    # plot
    plot<-ggplot(long_df, aes(key,y)) +
      geom_tile(aes(fill = z)) + 
      scale_fill_gradient2(low = "darkblue",
                           mid = "white",
                           high = "darkred",
                           midpoint = 0.99,
                           limits =c(0.98, 1))+
      #geom_text(aes(label = format(round(z, digits=3), nsmall = 3)),size=2.75) +
      ylab("") +  
      xlab("") + 
      theme_classic() + 
      theme(axis.text.x = element_text(size = 4, angle = 45, hjust = 1, color = "black"),
            axis.text.y = element_text(size = 4, color = "black"),
            axis.ticks.y = element_line(linewidth = 0.5),
            axis.ticks.x=element_line(linewidth = 0.5)) 
    
    # rows <- as.data.frame(wide_matrix)
    # height=nrow(rows) *0.15
    # width=height*1.2
    
    # save plot
    ggsave(file = paste0("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/15_DE_REPLICATION/Men_Gon_clusters/transformed_cluster_",clusterID,"_ANI_heatmap.pdf"), plot = plot, width = 11, height = 8, limitsize = FALSE)

}


```


The chosen N. gonorrhoeae genome for cluster 000002 shares a minimum of 98.91% ANI and a max of 99.85% with all other genomes in it's respective cluster. 


###### 1.11.4 Dereplication @ 99% ANI for N. gonorrhoaeae
```{bash, eval=FALSE}

mainDIR="/workspace/jmarkwelchlab"
projectID="P_0003_Neisseriaceae"
OUT=$mainDIR/$projectID/15_DE_REPLICATION/meningitidis_gonorrhoeae_ANI_99_dereplication_fastANI
DIR_ANI=$mainDIR/$projectID/15_DE_REPLICATION/meningitidis_gonorrhoeae_ANI_98_dereplication_fastANI/SIMILARITY_SCORES

anvi-dereplicate-genomes -o $OUT --ani-dir $DIR_ANI --skip-fasta-report --program fastANI --similarity-threshold 0.99 --cluster-method simple_greedy --representative-method centrality --num-threads 20 --log-file $mainDIR/$projectID/LOGS/meningitidis_gonorrhoeae_ANI_99_dereplication_fastANI.log

#Check
cat $OUT/CLUSTER_REPORT.txt | wc -l

# send to local
scp -r BARHAL_SERVER:/workspace/jmarkwelchlab/P_0003_Neisseriaceae/15_DE_REPLICATION/meningitidis_gonorrhoeae_ANI_99_dereplication_fastANI /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/15_DE_REPLICATION/
```



###### 1.11.5 Annotate with marker gene DBs (Bac_71, rRNAs Arch and Bac); tRNAs; GTDB taxonomy

```{bash, eval=FALSE update contigs dbs}
# Get the list of 21 genome IDs for the slected genomes
awk 'NR==FNR{a[$0]; next} $1 in a' DATA/Men_gon_selected_pangenome_IDs.txt $MEN_GNO_nameConversions | awk '{print $2}'> DATA/Men_gon_selected_genome_IDs.txt

for x in `cat DATA/Men_gon_selected_genome_IDs.txt`
do
cp 04_CONTIGS_DB/MEN_GNO_GENOMES/$x-contigs.db 04_CONTIGS_DB/MEN_GNO_GENOMES/$x-v20-contigs.db
anvi-migrate --migrate-safely 04_CONTIGS_DB/MEN_GNO_GENOMES/$x-contigs.db
done

# run hmms
for genome in `cat DATA/Men_gon_selected_genome_IDs.txt`
do
anvi-run-hmms -c 04_CONTIGS_DB/MEN_GNO_GENOMES/${genome}-contigs.db -T 10
done

```

We need to add these genomes to the master set, but the genome IDs need to be changed so that they don't overlap with the other non-men and non-gon genomes (i.e., the G_0615 set).

We can copy them into the main 04_CONTIGS_DB/ directory with their changed names.

First I used nano to manually edit the DATA/Men_gon_selected_genome_IDs.txt to include a column with new genome IDs ranging from 621 to 643. I saved the new file as DATA/Men_gon_selected_new_genome_IDs.txt. There is a header  old_ID and new_ID


```{bash, eval=FALSE}

# We need to change the project_names to reflect new genome IDs 
for x in `cat DATA/Men_gon_selected_genome_IDs.txt`
do
new_ID=$(cat DATA/Men_gon_selected_new_genome_IDs.txt | grep "$x" | awk '{print $2}')
cp 04_CONTIGS_DB/MEN_GNO_GENOMES/$x-contigs.db 04_CONTIGS_DB/$new_ID-contigs.db
anvi-db-info --self-key project_name --self-value "$new_ID" 04_CONTIGS_DB/$new_ID-contigs.db --just-do-it
done

```


We also need to update the fasta files with the new genome IDs so that we can build the bac 71 phylogeny
```{bash}

for x in `cat DATA/Men_gon_selected_genome_IDs.txt`
do
new_ID=$(cat DATA/Men_gon_selected_new_genome_IDs.txt | grep "$x" | awk '{print $2}')
cp 03_GENOMES_EDITED/MEN_GNO_GENOMES/$x.fa 03_GENOMES_EDITED/$new_ID.fa
cp 03_GENOMES_EDITED/MEN_GNO_GENOMES/$x.report.tsv 03_GENOMES_EDITED/$new_ID.report.tsv
done

```



###### 1.11.6 Add Men/Gon genomes to G_0615 dataset to make a G_0636 dataset


We need to add the genomes to the contig paths file:
```{bash, eval=FALSE }

# copy contigs path file
GENOMES_decontaminated=$mainDIR/$projectID/DATA/$projectID-decontaminated-contig_paths.txt
GENOMES_G_0636=$mainDIR/$projectID/DATA/$projectID-G_0636-contig_paths.txt
cp $GENOMES_decontaminated $GENOMES_G_0636

# Add new genomes to file
#1 . update name conversions file 

new_IDs=DATA/Men_gon_selected_new_genome_IDs.txt 
selected_pan_IDs=DATA/Men_gon_selected_pangenome_IDs.txt
MEN_GNO_nameConversions=DATA/meningitidis_gonorrhoeae-name_conversions.txt
selected_name_conv=DATA/meningitidis_gonorrhoeae-selected_name_conversions.txt
selected_new_name_conv=DATA/meningitidis_gonorrhoeae-selected_new_name_conversions.txt

# Get pangenome and genome IDs for selected genomes only.
awk 'NR==FNR{a[$0]; next} $1 in a' $selected_pan_IDs $MEN_GNO_nameConversions | awk '{print $1,$2}'> $selected_name_conv

# Change the second column of the first file "$selected_name_conv" so that the old_ID is replaced with the new_ID for each respective row.
awk 'NR==FNR {id[$1]=$2; next} $2 in id {$2=id[$2]} 1' OFS=' ' "$new_IDs" "$selected_name_conv" | awk 'BEGIN {OFS="\t"} {$1=$1} 1' > $selected_new_name_conv


# 2. make contigs path file
selected_new_name_conv=DATA/meningitidis_gonorrhoeae-selected_new_name_conversions.txt
new_IDs=DATA/Men_gon_selected_new_genome_IDs.txt 
MEN_GON_GENOMES=DATA/meningitidis_gonorrhoeae_selected-contig_paths.txt
echo -e 'name\tcontigs_db_path' > $MEN_GON_GENOMES
for genome in `cat $new_IDs | awk 'NR>1{print $2}'`
do
cat $selected_new_name_conv | grep "$genome" | awk -v workPath="$PWD" 'BEGIN{FS=OFS="\t"}{print $1,workPath"/04_CONTIGS_DB/"$2"-contigs.db"}' >> $MEN_GON_GENOMES
done

# 3. add new selected men/gon contigs path file to the G_0615 path file
GENOMES_decontaminated=DATA/${projectID}-decontaminated-contig_paths.txt
MEN_GON_GENOMES=DATA/meningitidis_gonorrhoeae_selected-contig_paths.txt
GENOMES_G_0636=DATA/${projectID}-G_0636-contig_paths.txt
cat $GENOMES_decontaminated <(tail -n +2 $MEN_GON_GENOMES) > $GENOMES_G_0636


# make ITEMS file for meningitidis_gonorrhoeae_selected
ITEMS_MEN_GON=DATA/meningitidis_gonorrhoeae_selected-add_info.items.txt
seqidHOMD=/workspace/jmarkwelchlab/HOMD_INFO/SEQID_info.txt
selected_new_name_conv=DATA/meningitidis_gonorrhoeae-selected_new_name_conversions.txt
echo -e "item\tGenome_ID\tSpecies\tAssembly_ID\tAssembly_ID_norm\tStrain\tGenome_in_HOMD\tHOMD_ID\tG_ID" > $ITEMS_MEN_GON
for gcaID in `cat $selected_new_name_conv | cut -f1 | awk -F'_id_' '{print $2}' | rev | sed -e 's/_/./' | rev `
do
gcaID_ed=$(echo $gcaID | sed -e 's/\./_/')
itemsID=$(grep "$gcaID_ed" $selected_new_name_conv | cut -f1 | awk -F'_str_' -v OFS="\t" '{print $0,$0,$1}')
strainID=$( grep "$gcaID_ed" $selected_new_name_conv | cut -f1 | awk -F'_str_' -v OFS="\t" '{print $2}' | awk -F'_id_' -v OFS="\t" '{print $1}' )
genomeHOMD=$( grep "$gcaID" $seqidHOMD)
if [ -z "$genomeHOMD" ]
then
inHOMD="no"
hmtID="NA"
else
inHOMD="yes"
hmtID=$( echo "$genomeHOMD"| awk 'BEGIN{FS=OFS="\t"}{print "HMT_"$2}')
fi
gID=$( grep "$gcaID_ed" $nameConversions | cut -f2 )
echo -e "$itemsID\t$gcaID_ed\t$gcaID\t$strainID\t$inHOMD\t$hmtID\t$gID" >> $ITEMS_MEN_GON
done

# concatenated men/gon items file with the G_0620 items file
ITEMS=${projectID}-add_info.items.txt
ITEMS_MEN_GON=DATA/meningitidis_gonorrhoeae_selected-add_info.items.txt
G_0636_ITEMS=DATA/G_0636-add_info.items.txt
cat $ITEMS <(tail -n +2 $ITEMS_MEN_GON) > $G_0636_ITEMS


# make new meta data file to add to pangenome
layersADD_G_0636=DATA/$projectID-G_0636-add_info.items.txt
G_0636_ITEMS=DATA/G_0636-add_info.items.txt
head -n1 $G_0636_ITEMS > $layersADD_G_0636
for genomeTypeID in `cat $GENOMES_G_0636`
do
cat $G_0636_ITEMS | grep "$genomeTypeID" >> $layersADD_G_0636
done

# Update name conversions
selected_new_name_conv=DATA/meningitidis_gonorrhoeae-selected_new_name_conversions.txt
G_0615_name_conv=DATA/G_0615-name_conversions.txt
G_0636_name_conv=DATA/G_0636-name_conversions.txt
# Add selected genomes listed above to name conversion file
cat $G_0615_name_conv $selected_new_name_conv > $G_0636_name_conv
```

## 1.12. GTDB classification

```{bash, eval=FALSE}


conda deactivate # ignore it if it complains
module purge
module load miniconda/3
source /bioware/miniconda3/bashrc 
conda activate /bioware/gtdbtk-2.3.0
gtdbtk check_install

module load clusters/barhal
module load jbpc


# Create the directory.

mkdir 29_GTDBTK/G_0636
mkdir 29_GTDBTK/G_0636/Genomes
mkdir 29_GTDBTK/G_0636/classify_out 
mkdir 29_GTDBTK/G_0636/tmp


projectID="P_0003_Neisseriaceae"
mainDIR="/workspace/jmarkwelchlab"

cd $mainDIR/$projectID

# make custom_taxonomy_file for gtdbtk decorate
cat DATA/$projectID-G_0636-add_info.items.txt | awk -F"\t" 'NR>1 {print $9, $2}' > 29_GTDBTK/G_0636/custom_taxonomy_file.txt


G_0636_name_conv=DATA/G_0636-name_conversions.txt
# move genomes to 29_GTDBTK/Genomes directory (moved to qrsh screen session to run this code)
for G_ID in `cat $G_0636_name_conv | awk '{print $2}'`; do
cp 03_GENOMES_EDITED/${G_ID}.fa 29_GTDBTK/G_0636/Genomes/${G_ID}.fa
done

# run script
# Added --mash_db gtdb-tk-r214.msh 
#chmod +x CRIPTS/script-gtdbtk-classify_wf_G_0636.sh

/workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/script-gtdbtk-classify_wf_G_0636.sh > LOGS/gtdbtk-classify_wf_G_0636.log 2>&1 &


```

```{r}
library("dplyr")
gtdb_results <- read.table("/workspace/jmarkwelchlab/P_0003_Neisseriaceae/29_GTDBTK/G_0636/classify_out/gtdbtk.bac120.summary.tsv", header = TRUE, sep = "\t")

meta <- read.table("/workspace/jmarkwelchlab/P_0003_Neisseriaceae/DATA/G_0636-name_conversions.txt", header = FALSE)
colnames(meta) <- c("pangenome_ID", "G_ID")

gtdb_results_merged <- merge(gtdb_results, meta, by.x = "user_genome", by.y = "G_ID")

gtdb_results_merged <- gtdb_results_merged %>% 
  mutate(gtdb_species_classification = classification)

# Remove all before and up to ";s_"
gtdb_results_merged$gtdb_species_classification <- gsub(".*;s__","",gtdb_results_merged$gtdb_species_classification)

# reorder columns
gtdb_results_merged <- gtdb_results_merged %>%
  select(user_genome,pangenome_ID,classification,gtdb_species_classification,classification_method, fastani_reference,fastani_reference_radius,fastani_taxonomy, fastani_ani,fastani_af, note, other_related_references.genome_id.species_name.radius.ANI.AF.)

# write results
write.table(gtdb_results_merged, "/workspace/jmarkwelchlab/P_0003_Neisseriaceae/29_GTDBTK/G_0636/classify_out/gtdbtk.bac120.merged_results.txt", quote = FALSE, row.names = FALSE, sep = "\t")
```


```{bash, eval=FALSE}

# send results to local
scp -r BARHAL_SERVER:/workspace/jmarkwelchlab/P_0003_Neisseriaceae/29_GTDBTK/G_0636/classify_out/gtdbtk.bac120.merged_results.txt /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/29_GTDBTK/G_0636/gtdbtk.bac120.merged_results.txt

```


# 2. Genome Dereplication

##### 2.1 G_0636

```{bash}

mainDIR="/workspace/jmarkwelchlab"
projectID="P_0003_Neisseriaceae"
DIR_Derep=$mainDIR/$projectID/15_DE_REPLICATION
DIR_ANI=$mainDIR/$projectID/09_PANGENOME/P_0003_Neisseriaceae_G_0636-pangenome/P_0003_Neisseriaceae_G_0636-pangenome-RESULTS/ANI_RESULTS

anvi-dereplicate-genomes -o $DIR_Derep/G_0636/G_0636_ANI_98_dereplication_fastANI --ani-dir $DIR_ANI --skip-fasta-report --program fastANI --similarity-threshold 0.98 --cluster-method simple_greedy --representative-method centrality --num-threads 20 --log-file $mainDIR/$projectID/LOGS/G_0636_ANI_98_dereplication_fastANI.log

#Check
cat $DIR_Derep/G_0636/G_0636_ANI_98_dereplication_fastANI/CLUSTER_REPORT.txt | wc -l

# send to local

scp -r BARHAL_SERVER:/workspace/jmarkwelchlab/P_0003_Neisseriaceae/15_DE_REPLICATION/G_0636/G_0636_ANI_98_dereplication_fastANI/CLUSTER_REPORT.txt /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/15_DE_REPLICATION/G_0636/CLUSTER_REPORT.txt

```

Beyond dereplication at 98% ANI, the follwoing criteria were also applied:

Remove 3 sp. genomes that associate with non-human-associated taxa.
2 with N. musculi
1 with N. zoodegmatis

Keep 3 N. sp. MVDL genomes that clustered in phylogeny next to N. anamolaris, but were distinct group based on ANI.

Remove the 3 N. macacae genomes and keep the 1 misidentified in NCBI as N. mucosa.

Remove 4 N. dentiae genomes and keep the 1 isolated from humans.

Remove N. zalophi  genome that is likely N. cinerea.

Keep the 2 N. perflava genomes that cluster in phylogeny far away from other N. perflava genomes because they are unique (share less than 90% ANI with all other genomes).

Keep both distinct K. sp. and N. sp. groups (2 genomes each).

Keep the K. potus genome that was isolated from humans - although likely an animal associated taxa. Thi sspecies was reported at very low rel abundance levels by metaphlan.Likely that the actual bacteria present is N. baccilliformis


##### 2.2 Plot distribution of selected reference genomes for mapping

```{r, eval=FALSE}

library("dplyr")

G_0213_derep_sel <- read.csv("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/DATA/Neisseriaceae_metapangenomics_mapping_reference_genomes.csv", header = TRUE)


G_0213_derep_sel_summary <- G_0213_derep_sel %>% 
  group_by(Genus,Reclassified_Species) %>% 
  summarise(n = n())


G_0213_derep_sel_summary <- as.data.frame(G_0213_derep_sel_summary)
G_0213_derep_sel_summary$Genus <- as.factor(G_0213_derep_sel_summary$Genus)
G_0213_derep_sel_summary$Reclassified_Species <- as.factor(G_0213_derep_sel_summary$Reclassified_Species)

dodge <- position_dodge(width=0.9)

G_0213_derep_sel_summary$Genus <- factor(G_0213_derep_sel_summary$Genus, levels = c("Neisseria", "Kingella", "Eikenella", "Simonsiella"), ordered = TRUE)


plot <- ggplot(G_0213_derep_sel_summary, aes(x = reorder(Reclassified_Species, -n), y = n, fill = Genus)) +
  geom_bar(stat = 'identity',position = dodge, color = "black") +
  geom_text(aes(label = n), vjust = -.5) +
  facet_grid(. ~ Genus, scales = "free", space = "free") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 45)) +
  xlab(NULL)+
  ylab("Frequency") +
  theme_classic() +
  theme(text = element_text(size = 12, color = "black"),
        strip.text.x = element_blank(),
        axis.text.x = element_text(angle = 45,vjust = 1, hjust=1, face = "italic", color = "black")) 

# clean up plot slightly
plot <- plot +
  theme(
    panel.background = element_rect(fill = "white"), 
    panel.grid = element_blank(),
    axis.line = element_blank()
  )



ggsave("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/FIGURES/Plots/Renamed_genomes_G_0213_dereplicated_selected.pdf", plot, height = 6, width = 11)




G_0213_derep_sel_summary <- G_0213_derep_sel %>% 
  group_by(Genus,Species) %>% 
  summarise(n = n())


G_0213_derep_sel_summary <- as.data.frame(G_0213_derep_sel_summary)
G_0213_derep_sel_summary$Genus <- as.factor(G_0213_derep_sel_summary$Genus)
G_0213_derep_sel_summary$pecies <- as.factor(G_0213_derep_sel_summary$Species)

dodge <- position_dodge(width=0.9)

G_0213_derep_sel_summary$Genus <- factor(G_0213_derep_sel_summary$Genus, levels = c("Neisseria", "Kingella", "Eikenella", "Simonsiella"), ordered = TRUE)


plot <- ggplot(G_0213_derep_sel_summary, aes(x = reorder(Species, -n), y = n, fill = Genus)) +
  geom_bar(stat = 'identity',position = dodge, color = "black") +
  geom_text(aes(label = n), vjust = -.5) +
  facet_grid(. ~ Genus, scales = "free", space = "free") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 45)) +
  xlab(NULL)+
  ylab("Frequency") +
  theme_classic() +
  theme(text = element_text(size = 12, color = "black"),
        strip.text.x = element_blank(),
        axis.text.x = element_text(angle = 45,vjust = 1, hjust=1, face = "italic", color = "black")) 

# clean up plot slightly
plot <- plot +
  theme(
    panel.background = element_rect(fill = "white"), 
    panel.grid = element_blank(),
    axis.line = element_blank()
  )



ggsave("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/FIGURES/Plots/NCBI_names_G_0213_dereplicated_selected.pdf", plot, height = 6, width = 11)
```


##### 2.3. Update reference genome metadata files

Send mapping reference genome metadata to Barhal
```{bash, eval=FALSE}

scp -r /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/DATA/Neisseriaceae_metapangenomics_mapping_reference_genomes.csv BARHAL_SERVER:/workspace/jmarkwelchlab/P_0003_Neisseriaceae/DATA/Neisseriaceae_metapangenomics_mapping_reference_genomes.csv
```

```{bash, eval=FALSE}

# Update contigs path file
GENOMES_G_0636=$mainDIR/$projectID/DATA/$projectID-G_0636-contig_paths.txt
GENOMES_G_0213=$mainDIR/$projectID/DATA/$projectID-G_0213-contig_paths.txt
mapping_ref_genome_metadata=$mainDIR/$projectID/DATA/Neisseriaceae_metapangenomics_mapping_reference_genomes.csv
mapping_ref_genomes=$mainDIR/$projectID/DATA/Pangenome_IDs_Neisseriaceae_metapangenomics_mapping_reference_genomes.txt

awk -F"," 'NR>1{print $1}' "$mapping_ref_genome_metadata" > "$mapping_ref_genomes"
grep -Ff "$mapping_ref_genomes" "$GENOMES_G_0636" > "$GENOMES_G_0213"


# Update meta data file to add to pangenome
layersADD_G_0213=$mainDIR/$projectID/DATA/$projectID-G_0213-add_info.items.txt
ITEMS=$mainDIR/$projectID/DATA/$projectID-G_0636-add_info.items.txt
head -n1 $ITEMS > $layersADD_G_0213
for genomeTypeID in `cat $GENOMES_G_0213`
do
cat $ITEMS | grep "$genomeTypeID" >> $layersADD_G_0213
done


# update name conversions file
G_0636_name_conv=$mainDIR/$projectID/DATA/G_0636-name_conversions.txt
G_0213_name_conv=$mainDIR/$projectID/DATA/G_0213-name_conversions.txt
grep -Ff "$mapping_ref_genomes" "$G_0636_name_conv" > "$G_0213_name_conv"

# make genomes ID file
G_0213_genomesID=$mainDIR/$projectID/DATA/id_genomes-G_0213.txt
awk '{print $2}' $G_0213_name_conv > $G_0213_genomesID

```

# 3. Concatenate selected reference genomes, make binning and decomposition files.


```{bash, eval=FALSE Concatenate fasta for reference genome set for mapping }

projectID="P_0003_Neisseriaceae"
mainDIR="/workspace/jmarkwelchlab"
DIR_Assemblies=02_ASSEMBLIES
DIR_Contigs=03_GENOMES_EDITED
DIR_ContigsDB=04_CONTIGS_DB
G_0213_genomesID=$mainDIR/$projectID/DATA/id_genomes-G_0213.txt
RAW_ASSEMBLY=$DIR_Assemblies/${projectID}-RAW.fa
PROJECT_CONTIGS=$DIR_Contigs/${projectID}.fa
PROJECT_REPORT=$DIR_Contigs/${projectID}.report.tsv
BINNING=$DIR_Contigs/${projectID}.binning.tsv
DECOMPOSE=$DIR_Contigs/${projectID}.decompose.tsv
CONTIGS_DB=$DIR_ContigsDB/${projectID}-contigs.db
minContigSIZE=300
G_0213_name_conv=$mainDIR/$projectID/DATA/G_0213-name_conversions.txt


# make concatenated genome file for METAPANGENOME
for genome in `cat $G_0213_genomesID`
do
CONTIGS=$DIR_Contigs/${genome}.fa
cat $CONTIGS >> $RAW_ASSEMBLY
done

# Reformat concatenated-genomes using ANVIO (anvi-script-reformat-fasta)
anvi-script-reformat-fasta -l $minContigSIZE -o $PROJECT_CONTIGS --simplify-names --prefix ${projectID} -r $PROJECT_REPORT $RAW_ASSEMBLY
```


To produce the decomposition file, please see Modifying_decomposition_file.Rmd. Below is an exmaple, but does not represent the process used to generate the modified decomposition file used in the decompose and summarize script.

```{bash, eval=FALSE}

# Binning and decomposition files
cat $PROJECT_REPORT | awk -F'_' -v OFS="_" 'NF{--NF};1' > $BINNING

while IFS= read -r line
do
intGenomeID=$( echo "$line" | awk '{print $2}' )
newBinID=$( echo "$line" | awk '{print $1}' )
grep "$intGenomeID" $BINNING | sed -e "s/$intGenomeID/$newBinID/" >> $DECOMPOSE
done < $G_0213_name_conv


```

# 4. Make concatenated Contig db, gene calling and hmms annotation

The following is a script for creating an Anvio contigs db (inlcuding gene calling and hmms annotation within anvio) for the concatenated reference genome set that we used for mapping. We need the single contigs db for the profiling steps; note that we have yet to annotate the contig db with gene functions, which we will do later. 

Note that running multiple Anvi'o processes that access and modify the same contigs database simultaneously can be risky, as concurrent writes or updates to the database may lead to corruption or data inconsistencies. This is especially true if one of the processes is modifying the database structure or adding new data to it, as is the case with anvi-run-ncbi-cogs. The anvi-run-ncbi-cogs function, as part of its operation, updates the contigs database with COG (Clusters of Orthologous Groups) functional annotations. This process involves writing new information to the database. On the other hand, anvi-profile generates profile databases from BAM files using a contigs database, and it might also access the database for read retrieval and other operations.

To avoid potential conflicts or database corruption, it's generally safer not to run these two commands simultaneously on the same contigs database. Here are a few suggestions:


```{bash}
# script for contigsdb, gene calling and hmms annotation (within anvio); need the single contigs db for the profiling steps; note that we have yet to annotate the contig db with gene functions. since this is just the first mapping run. 

projectID="P_0003_Neisseriaceae"

clusterize -n 10 -m jgiacomini@forsyth.org -log LOGS/Concatenated_contigdb.log /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/Mapping_profiling_scripts/script-01-contigsDB_annotation_hmms.sh $projectID 
```


# 5. Pangenome + ANI + Phylogeny

Pangenome: 
Using Anvi'o, we constructed a pangenome following previously developed methods (1, 3, 4, 18).  First, we used anvi-script-reformat-fasta to replace non-canonical nucleotide letters with "N" and remove from each reference genome all contigs with length less than 300 nt. We then converted each genome into an Anvi'o-compatible contigs database using anvi-gen-contigs-db. Open reading frames (ORFs), hereafter referred to as genes, were identified by Prodigal v2.6.3. Functional annotation of genes was achieved using multiple Anvi'o scripts, including anvi-run-hmms to find bacterial single-copy genes (Bacteria71 SCG set) (23, 24) with hidden Markov Model (HMM) profiles, anvi-run-ncbi-cogs using blastp (v2.10.1+) to annotate with the cluster of orthologous genes (COGs) database (version COG20) (25), and anvi-run-pfams and anvi-run-kegg-kofams with hmmscan from HMMER (v3.3.1) (26) to functionally annotate with Pfams (v34.0) (27) and KOfams/ KEGG Modules (v97.0) (28), respectively. We then used anvi-pan-genome to construct the annotated pangenome using BLASTP to calculate the amino acid-level identity between all possible gene pairs, with weak matches removed using the --minbit criterion of 0.5. The anvi-pan-genome program uses a Markov Cluster Algorithm (MCL) to group genes into putatively homologous groups called "gene clusters". We set the mcl-inflation parameter to 10 as suggested by Anvi’o for comparing very closely related genomes (https://merenlab.org/2016/11/08/pangenomics-v2/). Amino acid sequences within gene clusters were aligned with MUSCLE (v3.8.1551) (29). Finally, we performed hierarchical clustering across gene clusters and genomes using Euclidean distance and Ward linkage. This resulted in a pangenome showing the distribution of core and accessory genes across the reference genomes. 


We constructed a phylogenetic tree based on the amino acid sequences of 71 bacterial single-copy core genes (23, 24). We first used the Anvi’o program anvi-get-sequences-for-hmm-hits to align protein sequences using MUSCLE (v3.8.1551) (29), concatenate gene sequences, return only the most significant hit, and output amino acid sequences. Only genes that occurred in at least 50% of the genomes were used for the analysis, which in this case included all 71 genes. We trimmed alignments using trimAl (30) with the setting ‘-gt 0.5’ to remove all positions that were gaps in more than 50% of sequences. Maximum likelihood phylogenetic trees were then computed using IQ-TREE (31) with the WAG model (32) and 1000 bootstrap replicate support. We included a type strain genome for Burkholderia cepacia (strain BC16; GCA_009586235.1) to root the trees. To estimate pairwise whole genome average nucleotide identity (ANI) between the selected reference genomes in the pangenome we used the Anvi’o program anvi-compute-genome-similarity with the parameters ‘--program pyANI’ and ‘--method ANIb’. To compare genomes against the classification in GTDB, we used GTDB-Tk (version 2.3.0) (33) with classify_wf and the R214 reference data release. 


##### 5.1 Annotate contig data bases with functions

Pfam database version ........................: 36.0 (2023-07)

```{bash, eval=FALSE}

clusterize -n 20 -m jgiacomini@forsyth.org -log LOGS/G_0213_contigs_functional_annotations.log /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/Mapping_profiling_scripts/script-01-contigsDB_annotation_pfams_cogs.sh P_0003_Neisseriaceae

```


#####  5.2 Pangenome

```{bash, eval=FALSE}

export BLASTDB_LMDB_MAP_SIZE=10000000000

projectID="P_0003_Neisseriaceae"
mainDIR="/workspace/jmarkwelchlab"
COLLECTION=Genomes
site="SA"

DIR_MergedPROF=$mainDIR/$projectID/07_MERGED_PROFILE
DIR_DetectionGENOMES=$mainDIR/$projectID/13_DETECTED_GENOMES
DIR_SummaryPROF=$mainDIR/$projectID/08_PROFILE_SUMMARY
DIR_Pangenome=$mainDIR/$projectID/09_PANGENOME
DIR_ContigsDB=$mainDIR/$projectID/04_CONTIGS_DB

PROFILE_DB=$DIR_MergedPROF/${projectID}_${site}/PROFILE.db
CONTIGS_DB=$DIR_ContigsDB/${projectID}-contigs.db
INTERNAL=$DIR_Pangenome/${projectID}_G_0213-pangenome/internal_${projectID}.txt
genomes_G_0213=$mainDIR/$projectID/DATA/id_genomes-G_0213.txt
G_0213_bin_IDs=$mainDIR/$projectID/DATA/G_0213-name_conversions.txt
GENOMES_STORAGE=$DIR_Pangenome/${projectID}_G_0213-pangenome/${projectID}-GENOMES.db
PAN_RESULTS=$DIR_Pangenome/${projectID}_G_0213-pangenome/${projectID}-RESULTS
PAN_DB=$PAN_RESULTS/${projectID}-PAN.db
ANIb_RESULTS=$DIR_Pangenome/${projectID}_G_0213-pangenome/ANIb-RESULTS
ADD_INFO=$mainDIR/$projectID/DATA/$projectID-G_0213-add_info.items.txt
PAN_LAYERS=$DIR_Pangenome/${projectID}_G_0213-pangenome/${projectID}-add_info.layers.tsv


# Make pangenome directory
mkdir $DIR_Pangenome/${projectID}_G_0213-pangenome

# Make internal genome info file to be used for genome data base later
echo -e 'name\tbin_id\tcollection_id\tprofile_db_path\tcontigs_db_path' > $INTERNAL
cat $G_0213_bin_IDs | awk -v collection="$COLLECTION" -v profile_db="$PROFILE_DB" -v contigs_db="$CONTIGS_DB" 'BEGIN{FS=OFS="\t"}{print $1,$1,collection,profile_db,contigs_db}' >> $INTERNAL

# build annotated genome data base
clusterize -n 20 -m jgiacomini@forsyth.org -l LOGS/G_0213_genomes_storage.log anvi-gen-genomes-storage -i $INTERNAL --gene-caller prodigal -o $GENOMES_STORAGE

# build annotated pangenome
clusterize -n 20 -m jgiacomini@forsyth.org -l LOGS/G_0213_pangenome.log anvi-pan-genome -g $GENOMES_STORAGE --align-with muscle --use-ncbi-blast --minbit 0.5 --mcl-inflation 10 -n $projectID -o $PAN_RESULTS --num-threads 20 --enforce-hierarchical-clustering --I-know-this-is-not-a-good-idea


# create genome metadata $PAN_LAYERS 
G_0213_pangenome_IDs=$mainDIR/$projectID/DATA/G_0213_pangenome_IDs.txt
awk '{print $1}' $G_0213_name_conv > $G_0213_pangenome_IDs
grep -w -F -f $G_0213_pangenome_IDs -e item $ADD_INFO > $PAN_LAYERS

# import annotated reference genome metadata into pangenome
anvi-import-misc-data -p $PAN_DB -t layers $PAN_LAYERS --just-do-it
```


##### 5.3 Extract pangenome gene-cluster freq. tree

```{bash, eval=FALSE}
PanDir=$mainDIR/$projectID/09_PANGENOME/P_0003_Neisseriaceae_G_0213-pangenome

anvi-export-misc-data -p $PAN_DB --target-data-table layer_orders -o $PanDir/Layers_order
cat $PanDir/Layers_order | grep 'gene_cluster frequencies'| awk -F"\t" '{print $3}' >  $PanDir/gene_cluster_frequencies_newick

```

Send to local for detection plot:

```{bash, eval=FALSE}

scp -r BARHAL_SERVER:/workspace/jmarkwelchlab/P_0003_Neisseriaceae/09_PANGENOME/P_0003_Neisseriaceae_G_0213-pangenome/ /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/G_0213/

```



##### 5.4 ANI for G_0213 pangenome

```{bash, eval=FALSE}

# add ANI (ANIb method) to annotated pangenome
#module load blast+/2.10.1
clusterize -n 20 -m jgiacomini@forsyth.org -l LOGS/G_0213_ANIb.log anvi-compute-genome-similarity -i $INTERNAL -o $ANIb_RESULTS --pan-db $PAN_DB --program pyANI --method ANIb --num-threads 20

```


##### 5.4 BAC_71 SCG Phylogeny

Make sure you are in the anvio-8 environment:

screen -S qrsh_sesh
qrsh -pe allslots 5

module load anvio/8-conda
module load clusters/barhal
module load jbpc

Manually add outgroup to contigs paths file, then extract and concatenate aa sequences for bacteria_71 collection: 

```{bash, eval=FALSE}
projectID="P_0003_Neisseriaceae"
mainDIR="/workspace/jmarkwelchlab"
DIR_Phylo=$mainDIR/$projectID/14_PHYLOGENOMICS
GENOMES_G_0213=$mainDIR/$projectID/DATA/$projectID-G_0213-contig_paths.txt
GENOMES_G_0213_with_outgroup=$mainDIR/$projectID/DATA/$projectID-G_0213-contig_paths_with_outgroup.txt
num_genomes=213
min_num_genomes=$( expr $num_genomes / 2 )

# manually add out-group to contigs path file
echo -e 'name\tcontigs_db_path' >> $GENOMES_G_0213_with_outgroup
cat $GENOMES_G_0213_with_outgroup $GENOMES_G_0213 >> $GENOMES_G_0213_with_outgroup
echo -e 'Burkholderia_cepacia_str_BC16_id_GCA_009586235_1\t/workspace/jmarkwelchlab/P_0003_Neisseriaceae/14_PHYLOGENOMICS/Burkholderia_cepacia_str_BC16_id_GCA_009586235_1.db' >> $GENOMES_G_0213_with_outgroup

anvi-get-sequences-for-hmm-hits --e $GENOMES_G_0213_with_outgroup --hmm-source Bacteria_71 --min-num-bins-gene-occurs $min_num_genomes --get-aa-sequences --concatenate-genes --return-best-hit --align-with muscle -o $DIR_Phylo/Neisseriaceae_G_0213_bac71_sequences_with_outgroup

```

Build trees for each group with trimal/IQtree:

```{bash, eval=FALSE}

# trimal removes all positions in the alignment with gaps in 50% or more of the sequences
trimal -in $DIR_Phylo/Neisseriaceae_G_0213_bac71_sequences_with_outgroup -out $DIR_Phylo/Neisseriaceae_G_0213_bac71_sequences_with_outgroup.clean.fa -gt 0.50 

# build ML tree with WAG model and bootstrap supoprt
iqtree -s $DIR_Phylo/Neisseriaceae_G_0213_bac71_sequences_with_outgroup.clean.fa -nt AUTO -m WAG -bb 1000 -o "Burkholderia_cepacia_str_BC16_id_GCA_009586235_1"

```

Plot tree in R

```{bash, eval=FALSE}
scp -r BARHAL_SERVER:/workspace/jmarkwelchlab/P_0003_Neisseriaceae/14_PHYLOGENOMICS/Neisseriaceae_G_0213_bac71_sequences_with_outgroup.clean.fa.contree /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/14_PHYLOGENOMICS/Neisseriaceae_G_0213_bac71_sequences_with_outgroup.clean.fa.contree
```


Scale bar indicates the number of substitutions per site.
```{r, eval= FALSE}
library(dendextend)
library(ape)
library("dplyr")
library(phytools)
library(phylogram)

MLtree<-ape::read.tree("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/14_PHYLOGENOMICS/Neisseriaceae_G_0213_bac71_sequences_with_outgroup.clean.fa.contree")

pdf("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/14_PHYLOGENOMICS/G_0213_ML_tree_with_outgroup_bac71_IQTree.pdf",
    width = 10, height = 12)
plot(MLtree, cex = 0.25, edge.width = 0.5)
nodelabels(text = MLtree$node.label,node=2:MLtree$Nnode+Ntip(MLtree),frame="none",adj=c(1.2,-0.5), cex = 0.2, col = "red")
add.scale.bar()
dev.off() 

```


Chronos and dropped outgroup:

```{r, eval= FALSE}

library(dendextend)
library(ape)
library("dplyr")
library(phytools)
library(phylogram)

MLtree<-ape::read.tree("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/14_PHYLOGENOMICS/Neisseriaceae_G_0213_bac71_sequences_with_outgroup.clean.fa.contree")
MLtree_dropped_outgroup <- drop.tip(MLtree, "Burkholderia_cepacia_str_BC16_id_GCA_009586235_1")
dend_MLtree_dropped_outgroup <- ape::chronos(MLtree_dropped_outgroup)

pdf("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/14_PHYLOGENOMICS/G_0213_ML_tree_no_outgroup_bac71_IQTree.pdf",
    width = 10, height = 12)
plot(dend_MLtree_dropped_outgroup, cex = 0.25, edge.width = 0.5)
nodelabels(text = MLtree_dropped_outgroup$node.label,node=2:MLtree_dropped_outgroup$Nnode+Ntip(MLtree_dropped_outgroup),frame="none",adj=c(1.2,-0.5), cex = 0.2, col = "red")
add.scale.bar()
dev.off() 


```

Save un-rooted tree as newick file that we can import into the pangenome:

```{r, eval=FALSE}

write.tree(dend_MLtree_dropped_outgroup, file = "/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/14_PHYLOGENOMICS/Neisseriaceae_G_0213_bac71_sequences_with_outgroup.clean.fa.contree.newick")

```

Send to Barhal:

```{bash, eval=FALSE}

scp -r /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/14_PHYLOGENOMICS/Neisseriaceae_G_0213_bac71_sequences_with_outgroup.clean.fa.contree.newick BARHAL_SERVER:/workspace/jmarkwelchlab/P_0003_Neisseriaceae/14_PHYLOGENOMICS/Neisseriaceae_G_0213_bac71_sequences_with_outgroup.clean.fa.contree.newick

```

##### 5.5 Add Bacteria 71 SCG phylogeny to pangenome

```{bash, eval=FALSE}

# Note that we built the SCG tree earlier
# import tree to pangenome
projectID="P_0003_Neisseriaceae"
mainDIR="/workspace/jmarkwelchlab"
DIR_Phylo=$mainDIR/$projectID/14_PHYLOGENOMICS
phyName_bac71=Bac71_phylogeny
phyloTree_bac71=$DIR_Phylo/Neisseriaceae_G_0213_bac71_sequences_with_outgroup.clean.fa.contree.newick
ADD_phyloTree_bac71=$DIR_Phylo/${phyName_bac71}.layer_orders.tsv
DIR_Pangenome=$mainDIR/$projectID/09_PANGENOME
PAN_RESULTS=$DIR_Pangenome/${projectID}_G_0213-pangenome/${projectID}-RESULTS
PAN_DB=$PAN_RESULTS/${projectID}-PAN.db

awk -v phyCount="$phyName_bac71" 'BEGIN{FS=OFS="\t"}NR==1{print "item_name","data_type","data_value"}{print phyCount,"newick",$0}' $phyloTree_bac71 > $ADD_phyloTree_bac71

anvi-import-misc-data -t layer_orders -p $PAN_DB $ADD_phyloTree_bac71

```

##### 5.4 b. PAN SCG Phylogeny

screen -S qrsh_sesh
qrsh -pe allslots 10

module load anvio
module load clusters/barhal
module load jbpc

```{bash, eval=FALSE}

projectID="P_0003_Neisseriaceae"
mainDIR="/workspace/jmarkwelchlab"
DIR_Pangenome=$mainDIR/$projectID/09_PANGENOME
PAN_DB=$DIR_Pangenome/${projectID}_G_0213-pangenome/${projectID}-RESULTS/${projectID}-PAN.db
GENOMES_DB=$DIR_Pangenome/${projectID}_G_0213-pangenome/${projectID}-GENOMES.db
N_genomes=213


anvi-get-sequences-for-gene-clusters -g $GENOMES_DB \
                                     -p $PAN_DB \
                                     -o $DIR_Pangenome/single_copy_core_genes-fasta \
                                     --max-num-genes-from-each-genome 1 \
                                     --min-num-genomes-gene-cluster-occurs $N_genomes \
                                     --concatenate-gene-clusters 
```

Anvio reported "Your filters resulted in 222 gene clusters that contain a total of 47286 genes.". This means that a total of 222 single copy core genes were shared among the 213 genomes.  

Build trees for each group with trimal/IQtree:
  
```{bash, eval=FALSE}

# trimal removes all positions in the alignment with gaps in 50% or more of the sequences
trimal -in $DIR_Pangenome/single_copy_core_genes-fasta -out $DIR_Pangenome/single_copy_core_genes-fasta.clean.fa -gt 0.50 

# build ML tree with WAG model and bootstrap supoprt
iqtree -s $DIR_Pangenome/single_copy_core_genes-fasta.clean.fa -nt AUTO -m WAG -bb 1000 

```

Send to local...
```{bash, eval=FALSE}
scp -r BARHAL_SERVER:/workspace/jmarkwelchlab/P_0003_Neisseriaceae/09_PANGENOME/single_copy_core_genes-fasta.clean.fa.contree /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/14_PHYLOGENOMICS/G_0213_Pangenome_single_copy_core_genes-fasta.clean.fa.contree
```


##### 5.6 Display pangenome

Send pangenome to local (easier to display)

```{bash, eval=FALSE}

scp -r BARHAL_SERVER:/workspace/jmarkwelchlab/P_0003_Neisseriaceae/09_PANGENOME/P_0003_Neisseriaceae_G_0213-pangenome/ /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/09_PANGENOME/G_0213/

```

Run on local!!

```{bash, eval=FALSE}

pan_dir=/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/09_PANGENOME/G_0213

anvi-display-pan -p $pan_dir/P_0003_Neisseriaceae-RESULTS/P_0003_Neisseriaceae-PAN.DB -g $pan_dir/P_0003_Neisseriaceae-GENOMES.db


```


I created a spaced and not_spaced version. Both are colored according to the G_0636 ANI figure.

##### 5.7 Re-extract custom pangeome gene-cluster freq tree

Extract from local pangenome with saved state files.

```{bash, eval=FALSE}
PanDir=/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/09_PANGENOME/G_0213

anvi-export-misc-data -p $PanDir/P_0003_Neisseriaceae-RESULTS/P_0003_Neisseriaceae-PAN.DB --target-data-table layer_orders -o $PanDir/Layers_order

# The following code will print the names of the different trees
cat $PanDir/Layers_order | awk -F"\t" '{print $1}'

# This code will extract the tree of choice into a new file
cat $PanDir/Layers_order | grep 'gene cluster freq custom'| awk -F"\t" '{print $3}' >  $PanDir/gene_cluster_frequencies_custom_newick

```


At this point I re-created the detection and coverage figures using the new genomes order file 'gene_cluster_frequencies_custom_newick'.

##### 5.8 Add gene cluster collection (Core, Acc. Singles)

We then generated gene cluster collections for Core: gene clusters shared between all genomes, Accessory: gene clusters shared between two or more genomes, but not all, and Singletons: gene clusters found only in a single genome.
```{bash make pangenome collection}
iDir=$DIR_Pangenome/internal_annotated_${projectID}
genomeDB=$iDir/${projectID}-GENOMES.db
panDB=$iDir/${projectID}-RESULTS/${projectID}-PAN.db
numGenomes="77"

# get all gene-clusters
anvi-get-sequences-for-gene-clusters -g $genomeDB -p $panDB --min-num-genomes-gene-cluster-occurs 1 -o $iDir/${projectID}-accessoryGC.temp.txt
# get genusCore gene-clusters
anvi-get-sequences-for-gene-clusters -g $genomeDB -p $panDB --min-num-genomes-gene-cluster-occurs $numGenomes -o $iDir/${projectID}-genusCoreGC.temp.txt
# get singleton gene-clusters
anvi-get-sequences-for-gene-clusters -g $genomeDB -p $panDB --max-num-genomes-gene-cluster-occurs 1 -o $iDir/${projectID}-singletonsGC.temp.txt

# make list with unique gene-cluster IDs for genusCore & singletons
# add bin name in second column
# append to cas3 collection file
cat $iDir/${projectID}-genusCoreGC.temp.txt | grep ">" | awk -F'|' '{print $2}' | sed -e 's/gene_cluster\://' | sort | uniq | awk 'BEGIN{FS=OFS="\t"}{print $1, "Genus_core"}' > $iDir/${projectID}-collection.cas3.txt
cat $iDir/${projectID}-singletonsGC.temp.txt | grep ">" | awk -F'|' '{print $2}' | sed -e 's/gene_cluster\://' | sort | uniq | awk 'BEGIN{FS=OFS="\t"}{print $1, "Singletons"}' >> $iDir/${projectID}-collection.cas3.txt

# remove genus core & singletons gene-cluster IDs from accessory file
# print only missing gene clusters from accessoryID
# add bin name in second column
# append to cas3 collection file
comm -23 <( cat $iDir/${projectID}-accessoryGC.temp.txt | grep ">" | awk -F'|' '{print $2}' | sed -e 's/gene_cluster\://' | sort | uniq ) <( cat $iDir/${projectID}-collection.cas3.txt | cut -f1 | sort | uniq) | awk 'BEGIN{FS=OFS="\t"}{print $1, "Accessory"}' >> $iDir/${projectID}-collection.cas3.txt

# make collection info file
echo -e "Genus_core\tUNKNOWN_SOURCE\t#870000" > $iDir/${projectID}-collection.cas3-info.txt
echo -e "Accessory\tUNKNOWN_SOURCE\t#f2f2f2" >> $iDir/${projectID}-collection.cas3-info.txt
echo -e "Singletons\tUNKNOWN_SOURCE\t#db9e04" >> $iDir/${projectID}-collection.cas3-info.txt

# import cas3 collection & info to pangenome
anvi-import-collection -p $panDB -C cas3 --bins-info $iDir/${projectID}-collection.cas3-info.txt $iDir/${projectID}-collection.cas3.txt

rm $iDir/*temp.txt
```


# 6. FIGURE: ANI heatmap

Send results to local

```{bash, eval = FALSE}


scp -r BARHAL_SERVER:/workspace/jmarkwelchlab/P_0003_Neisseriaceae/09_PANGENOME/P_0003_Neisseriaceae_G_0213-pangenome/ANIb-RESULTS /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/ANI/G_0213

```



##### 6.1 Percent identity plot

```{r, eval=FALSE}

# packages
library("ggplot2")
library("tidyr")
library("reshape2")
library("phytools")
library("tibble")
library("dendextend")
library("ape")
library("dplyr")
library("ggtree")
library("cowplot")

# load ANI data
df <-read.table("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/ANI/G_0213/ANIb_percentage_identity.txt", header = TRUE)

# load ANI tree
df_newick <-  ape::read.tree("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/ANI/G_0213/ANIb_percentage_identity.newick")
p <- ggtree(df_newick) + geom_treescale() + geom_tiplab(align=TRUE, linetype = "dotted", linesize = 0.3, size = 0.5, colour = "black")
p 
ggsave("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/FIGURES/ANI/Test_tree.pdf", p, height = 10, width = 8)

  
#   orders <- phytools::compute.mr(df_newick, type = "matrix")
# df_newick_rownames <- rev(rownames(df_newick_orders))
# data_ordered <- df[ order(match(df$key, df_newick_rownames)), ]
# df_newick_rownames_df <- as.data.frame(df_newick_rownames)
# cols_A<- ncol(df) -1
# df_newick_rownames_df2 <- as.data.frame(matrix(df_newick_rownames_df$df_newick_rownames, ncol = cols_A, byrow = TRUE))
# names(df_newick_rownames_df2) <- df_newick_rownames_df2[1,]
# df_newick_rownames_df3 <- df_newick_rownames_df2[-1,]
# df_newick_rownames_df3 <- df_newick_rownames_df3 %>% 
#   add_column(key = NA, .before = 1)
# data_ordered2<-data_ordered[names(df_newick_rownames_df3)]

cols<- ncol(df) 
long_df <- reshape2::melt(data_ordered2,
                          id.vars=c("key"),
                          measure.vars=colnames(df[2:cols]),
                          variable.name="y",
                          value.name="z")
mylevels1 <- rev(get_taxa_name(p))
long_df$key <- factor(long_df$key,levels=mylevels1)
long_df$y <- factor(long_df$y, levels=mylevels1)



# get tick colors

# import colors df
# load csv file with genome IDs and colors; theis file was made when the tanglegram was made
bins_and_colors <-read.csv("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/13_DETECTED_GENOMES/pangenome_custom_gene_freq_order.csv", header = TRUE)

# Ensure 'Pangenome_ID' and 'key' are factors
bins_and_colors$Pangenome_ID <- factor(bins_and_colors$Pangenome_ID, levels = mylevels1)

# Sort the data frame by 'Pangenome_ID' using the defined factor levels
bins_and_colors <- bins_and_colors[order(bins_and_colors$Pangenome_ID), ]

# Extract 'Group_color' in the correct order
group_tick_colors <- rev(bins_and_colors$Group_color)

# plot
plot<-ggplot(long_df, aes(key,y)) +
  geom_tile(aes(fill = z)) + 
  scale_fill_gradient2(low = "darkblue",
                       mid = "white",
                       high = "darkred",
                       midpoint = 0.95,
                       limits =c(0.90, 1),
                       na.value="darkblue",
                        name = "ANI (%)")+
  #geom_text(aes(label = format(round(z, digits=3), nsmall = 3)),size=2.75) +
  ylab("") +  
  xlab("") + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 2),
        axis.text.y = element_text(color = "black", size = 2),
        axis.ticks.y = element_line(color = rev(group_tick_colors), linewidth = 1),
        axis.ticks.x = element_line(color = rev(group_tick_colors), linewidth = 1.3),
        axis.ticks.length=unit(.5, "cm")) 

# rows <- as.data.frame(df_newick$tip.label)
# height=nrow(rows) *0.3
# width=height*1.3
  
ggsave(file = "/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/ANI/G_0213/ANI_heatmap.pdf", plot = plot, width = 11, height = 8, limitsize = FALSE)


# Combine tree and heatmap plots
plot2 <- cowplot::plot_grid(p, plot, align = "h", rel_widths = c(1, 2))


ggsave(file = "/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/ANI/G_0213/ANI_heatmap_2.pdf", plot = plot2, width = 11, height = 8, limitsize = FALSE)


```

The distance bar represents the degree of genomic dissimilarity between sequences, as derived from pairwise Average Nucleotide Identity (ANI) values. Longer branches indicate greater dissimilarity, while shorter branches suggest closer similarity between genomes.


##### 6.2 Full percent identity plot

```{r, eval=FALSE}

# packages
library(ggplot2)
library(tidyr)
library(reshape2)
library(phytools)
library(tibble)
library(dendextend)
library(ape)
library("dplyr")

# load ANI data
df <-read.table("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/ANI/G_0213/ANIb_full_percentage_identity.txt", header = TRUE)

# load ANI tree
df_newick <-  ape::read.tree("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/ANI/G_0213/ANIb_full_percentage_identity.newick")

# order ANI data by ANI tree
df_newick_orders <- phytools::compute.mr(df_newick, type = "matrix")
df_newick_rownames <- rev(rownames(df_newick_orders))
data_ordered <- df[ order(match(df$key, df_newick_rownames)), ]
df_newick_rownames_df <- as.data.frame(df_newick_rownames)
cols_A<- ncol(df) -1
df_newick_rownames_df2 <- as.data.frame(matrix(df_newick_rownames_df$df_newick_rownames, ncol = cols_A, byrow = TRUE))
names(df_newick_rownames_df2) <- df_newick_rownames_df2[1,]
df_newick_rownames_df3 <- df_newick_rownames_df2[-1,]
df_newick_rownames_df3 <- df_newick_rownames_df3 %>% 
  add_column(key = NA, .before = 1)
data_ordered2<-data_ordered[names(df_newick_rownames_df3)]

cols2<- ncol(data_ordered2) 
long_df <- reshape2::melt(data_ordered2,
                          id.vars=c("key"),
                          measure.vars=colnames(data_ordered2[2:cols2]),
                          variable.name="y",
                          value.name="z")
mylevels1 <- df_newick$tip.label
long_df$key <- factor(long_df$key,levels=mylevels1)
long_df$y <- factor(long_df$y, levels=mylevels1)



# get tick colors

# import colors df
# load csv file with genome IDs and colors; theis file was made when the tanglegram was made
bins_and_colors <-read.csv("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/13_DETECTED_GENOMES/pangenome_custom_gene_freq_order.csv", header = TRUE)

# Ensure 'Pangenome_ID' and 'key' are factors
bins_and_colors$Pangenome_ID <- as.factor(bins_and_colors$Pangenome_ID)

# Match the order of 'Pangenome_ID' in bins_and_colors to the order in 'key' of long_df
order_indices <- match(long_df$key, bins_and_colors$Pangenome_ID)

# Use the order_indices to reorder bins_and_colors
bins_and_colors_reordered <- bins_and_colors[order_indices, ]


group_tick_colors <- bins_and_colors_reordered$Group_color


# plot
plot<-ggplot(long_df, aes(key,y)) +
  geom_tile(aes(fill = z)) + 
  scale_fill_gradient2(low = "darkblue",
                       mid = "white",
                       high = "darkred",
                       midpoint = 0.75,
                       limits =c(0.5, 1),
                       na.value="darkblue",
                        name = "ANI (%)")+
  #geom_text(aes(label = format(round(z, digits=3), nsmall = 3)),size=2.75) +
  ylab("") +  
  xlab("") + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 2),
        axis.text.y = element_text(color = "black", size = 2),
        axis.ticks.y = element_line(color = rev(group_tick_colors), size = 1),
        axis.ticks.x = element_line(color = rev(group_tick_colors), size = 1.3),
        axis.ticks.length=unit(.5, "cm")) 

# rows <- as.data.frame(df_newick$tip.label)
# height=nrow(rows) *0.3
# width=height*1.3
  
ggsave(file = "/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/ANI/G_0213/ANIb_full_percentage_identity_heatmap.pdf", plot = plot, width = 11, height = 8, limitsize = FALSE)

```

##### 6.3 Coverage 


Alignment coverage is the proportion of the query genome that aligns against the reference genome. This can be asymmetrical even when the alignment itself is symmetrical, as the genomes participating in a pairwise alignment may have differing amounts of genomic sequence that do not contribute to the alignment.

Percentage coverage matrix for Neisseria, Kingella, Eikenlla and Simonsiella ANIb analysis -  coverage is calculated as alignment length / genome length.

Each cell represents a pairwise comparison between the named genomes on rows and columns, and the number in each cell is the pairwise coverage of each genome by aligned regions in the comparison (i.e, Percentage coverage). The color scheme colors cells as white with a percentage coverage equal to 0.50, greater as red, and less than as s blue. This division corresponds to a strict majority of each genome in the comparison being align-able (a plausible ad hoc minimum requirement for two sequences being considered “the same thing”).


In our case, nearly all the pairwise comparisons between Eikenella species and Neisseria, Kingella and Simonsiella are blue cells below 50%. This indicates that the proportion of each genome that aligns is very small, and we are safe to assert that these organisms come from different genera. The exception to this is the comparison between BVAF and floridianus: their coverage is higher, at ≈15%. This may indicate a common plasmid or mobile element, or it may indicate a more recent common ancestor than the other comparisons; they may or may not validly be in the same genus - we would need to investigate further to understand their relationship.

The BPEN/640 comparison is conclusive, however. Their coverage/AF is essentially 100%, so these are closely-related, highly sequence-homologous organisms.

```{r, eval=FALSE}

# packages
library(ggplot2)
library(tidyr)
library(reshape2)
library(phytools)
library(tibble)
library(dendextend)
library(ape)
library("dplyr")

# load ANI data
df <-read.table("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/ANI/G_0213/ANIb_alignment_coverage.txt", header = TRUE)

# load ANI tree
df_newick <-  ape::read.tree("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/ANI/G_0213/ANIb_alignment_coverage.newick")

# order ANI data by ANI tree
df_newick_orders <- phytools::compute.mr(df_newick, type = "matrix")
df_newick_rownames <- rev(rownames(df_newick_orders))
data_ordered <- df[ order(match(df$key, df_newick_rownames)), ]
df_newick_rownames_df <- as.data.frame(df_newick_rownames)
cols_A<- ncol(df) -1
df_newick_rownames_df2 <- as.data.frame(matrix(df_newick_rownames_df$df_newick_rownames, ncol = cols_A, byrow = TRUE))
names(df_newick_rownames_df2) <- df_newick_rownames_df2[1,]
df_newick_rownames_df3 <- df_newick_rownames_df2[-1,]
df_newick_rownames_df3 <- df_newick_rownames_df3 %>% 
  add_column(key = NA, .before = 1)
data_ordered2<-data_ordered[names(df_newick_rownames_df3)]

cols2<- ncol(data_ordered2) 
long_df <- reshape2::melt(data_ordered2,
                          id.vars=c("key"),
                          measure.vars=colnames(data_ordered2[2:cols2]),
                          variable.name="y",
                          value.name="z")
mylevels1 <- df_newick$tip.label
long_df$key <- factor(long_df$key,levels=mylevels1)
long_df$y <- factor(long_df$y, levels=mylevels1)



# get tick colors

# import colors df
# load csv file with genome IDs and colors; theis file was made when the tanglegram was made
bins_and_colors <-read.csv("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/13_DETECTED_GENOMES/pangenome_custom_gene_freq_order.csv", header = TRUE)

# Ensure 'Pangenome_ID' and 'key' are factors
bins_and_colors$Pangenome_ID <- as.factor(bins_and_colors$Pangenome_ID)

# Match the order of 'Pangenome_ID' in bins_and_colors to the order in 'key' of long_df
order_indices <- match(long_df$key, bins_and_colors$Pangenome_ID)

# Use the order_indices to reorder bins_and_colors
bins_and_colors_reordered <- bins_and_colors[order_indices, ]


group_tick_colors <- bins_and_colors_reordered$Group_color


# plot
plot<-ggplot(long_df, aes(key,y)) +
  geom_tile(aes(fill = z)) + 
  scale_fill_gradient2(low = "darkblue",
                       mid = "white",
                       high = "darkred",
                       midpoint = 0.5,
                       limits =c(0, 1),
                       na.value="darkblue",
                       name = expression(atop("ANIb alignment coverage", "(Coverage/Aligned Fraction)"))) + 
  #geom_text(aes(label = format(round(z, digits=3), nsmall = 3)),size=2.75) +
  ylab("") +  
  xlab("") + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 2),
        axis.text.y = element_text(color = "black", size = 2),
        axis.ticks.y = element_line(color = rev(group_tick_colors), size = 1),
        axis.ticks.x = element_line(color = rev(group_tick_colors), size = 1.3),
        axis.ticks.length=unit(.5, "cm")) 

# rows <- as.data.frame(df_newick$tip.label)
# height=nrow(rows) *0.3
# width=height*1.3
  
ggsave(file = "/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/ANI/G_0213/ANIb_alignment_coverage_heatmap.pdf", plot = plot, width = 11, height = 8, limitsize = FALSE)

```

##### 6.4 Haddard matrix (Identity * coverage)

```{r, eval=FALSE}

# packages
library(ggplot2)
library(tidyr)
library(reshape2)
library(phytools)
library(tibble)
library(dendextend)
library(ape)
library("dplyr")

# load ANI data
df <-read.table("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/ANI/G_0213/ANIb_hadamard.txt", header = TRUE)

# load ANI tree
df_newick <-  ape::read.tree("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/ANI/G_0213/ANIb_hadamard.newick")

# order ANI data by ANI tree
df_newick_orders <- phytools::compute.mr(df_newick, type = "matrix")
df_newick_rownames <- rev(rownames(df_newick_orders))
data_ordered <- df[ order(match(df$key, df_newick_rownames)), ]
df_newick_rownames_df <- as.data.frame(df_newick_rownames)
cols_A<- ncol(df) -1
df_newick_rownames_df2 <- as.data.frame(matrix(df_newick_rownames_df$df_newick_rownames, ncol = cols_A, byrow = TRUE))
names(df_newick_rownames_df2) <- df_newick_rownames_df2[1,]
df_newick_rownames_df3 <- df_newick_rownames_df2[-1,]
df_newick_rownames_df3 <- df_newick_rownames_df3 %>% 
  add_column(key = NA, .before = 1)
data_ordered2<-data_ordered[names(df_newick_rownames_df3)]

cols2<- ncol(data_ordered2) 
long_df <- reshape2::melt(data_ordered2,
                          id.vars=c("key"),
                          measure.vars=colnames(data_ordered2[2:cols2]),
                          variable.name="y",
                          value.name="z")
mylevels1 <- df_newick$tip.label
long_df$key <- factor(long_df$key,levels=mylevels1)
long_df$y <- factor(long_df$y, levels=mylevels1)



# get tick colors

# import colors df
# load csv file with genome IDs and colors; theis file was made when the tanglegram was made
bins_and_colors <-read.csv("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/13_DETECTED_GENOMES/pangenome_custom_gene_freq_order.csv", header = TRUE)

# Ensure 'Pangenome_ID' and 'key' are factors
bins_and_colors$Pangenome_ID <- as.factor(bins_and_colors$Pangenome_ID)

# Match the order of 'Pangenome_ID' in bins_and_colors to the order in 'key' of long_df
order_indices <- match(long_df$key, bins_and_colors$Pangenome_ID)

# Use the order_indices to reorder bins_and_colors
bins_and_colors_reordered <- bins_and_colors[order_indices, ]


group_tick_colors <- bins_and_colors_reordered$Group_color


# plot
plot<-ggplot(long_df, aes(key,y)) +
  geom_tile(aes(fill = z)) + 
  scale_fill_gradient2(low = "darkblue",
                       mid = "white",
                       high = "darkred",
                       midpoint = 0.75,
                       limits =c(0.5, 1),
                       na.value="darkblue",
                       name = expression(atop("ANIb Hadamard", "(% Identity * Coverage)"))) +
  #geom_text(aes(label = format(round(z, digits=3), nsmall = 3)),size=2.75) +
  ylab("") +  
  xlab("") + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 2),
        axis.text.y = element_text(color = "black", size = 2),
        axis.ticks.y = element_line(color = rev(group_tick_colors), size = 1),
        axis.ticks.x = element_line(color = rev(group_tick_colors), size = 1.3),
        axis.ticks.length=unit(.5, "cm")) 

# rows <- as.data.frame(df_newick$tip.label)
# height=nrow(rows) *0.3
# width=height*1.3
  
ggsave(file = "/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/ANI/G_0213/ANIb_hadamard_heatmap.pdf", plot = plot, width = 11, height = 8, limitsize = FALSE)

```

# 7. FIGURE: Phylogeny and Pangenome tanglegram


##### 7.1. Tanglegram figure

```{r, eval=FALSE}

# load libraries
library(dendextend)
library(ape)
library("dplyr")
library(phytools)
library(phylogram)
library(gplots)

# load trees
PANtree<-ape::read.tree("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/09_PANGENOME/G_0213/gene_cluster_frequencies_custom_newick")

Bac_71_MLtree<-ape::read.tree("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/14_PHYLOGENOMICS/Neisseriaceae_G_0213_bac71_sequences_with_outgroup.clean.fa.contree")

# remove outgroup from tree file for merging into pangenome and making tanglegram
Bac_71_MLtree <- drop.tip(Bac_71_MLtree, "Burkholderia_cepacia_str_BC16_id_GCA_009586235_1")

# convert to dendrogram objects
dend_Bac_71_MLtree <- ape::chronos(Bac_71_MLtree)
dend_Bac_71_MLtree_2 <- as.dendrogram.phylo(dend_Bac_71_MLtree)
dend_PANtree <- ape::chronos(PANtree)
dend_PANtree_2 <- rev(as.dendrogram.phylo(dend_PANtree))
# load colors df
bins_and_colors <-read.csv("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/13_DETECTED_GENOMES/pangenome_custom_gene_freq_order.csv", header = TRUE)

# Ensure 'Pangenome_ID' and 'key' are factors
bins_and_colors$Pangenome_ID <- as.factor(bins_and_colors$Pangenome_ID)

#### Pangenome colors ####
# # Extract labels from pangenome dendrogram 
pan_labels <- dend_PANtree_2 %>% set("labels_to_char") %>% labels 

# Match the order of 'Pangenome_ID' in bins_and_colors to the order in 'key' of long_df
Pan_order_indices <- match(pan_labels, bins_and_colors$Pangenome_ID)

# Use the order_indices to reorder bins_and_colors
Pan_bins_and_colors_reordered <- bins_and_colors[Pan_order_indices, ]

# set color variable for plotting
Pan_group_colors <- Pan_bins_and_colors_reordered$Group_color


#### Phylo colors ####
# # Extract labels from bac71 phylo dendrogram 
bac71_phylo_labels <- dend_Bac_71_MLtree_2 %>% set("labels_to_char") %>% labels 

# Match the order of 'Pangenome_ID' in bins_and_colors to the order in 'key' of long_df
Phylo_order_indices <- match(bac71_phylo_labels, bins_and_colors$Pangenome_ID)

# Use the order_indices to reorder bins_and_colors
Phylo_bins_and_colors_reordered <- bins_and_colors[Phylo_order_indices, ]

# set color variable for plotting
Phylo_group_colors <- Phylo_bins_and_colors_reordered$Group_color

# make dendrogram list; second dedrogram will be fixed
bac71_SCG_PAN <- dendlist(dend_Bac_71_MLtree_2 %>% 
    set("labels_col", value = Phylo_group_colors),
    dend_PANtree_2 %>% 
      set("labels_col", value = Pan_group_colors)) 


# set height for pdf
rows <- as.data.frame(Bac_71_MLtree$tip.label)
height=nrow(rows) *0.0825

# build tanglegram
bac71_SCG_PAN_TANGLEGRAM <- bac71_SCG_PAN %>% dendextend::untangle(method="step1side") 

# set colors of lines; needs to be based on lefthand side dendrogram
bac71_PAN_TANGLEGRAM_labels <- bac71_SCG_PAN_TANGLEGRAM[[1]] %>%  labels
bac71_PAN_TANGLEGRAM_labels <- as.data.frame(bac71_PAN_TANGLEGRAM_labels)

line_colors <- merge(bac71_PAN_TANGLEGRAM_labels, Phylo_bins_and_colors_reordered, by.x="bac71_PAN_TANGLEGRAM_labels", by.y="Pangenome_ID", sort=F)
line_colors2 <- as.character(line_colors$Group_color)
line_colors3 <- col2hex(line_colors2)

# plot tanglegram
pdf(file = "/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/FIGURES/Phylogenies/G_0213_Phylogeny_vs_pangenome_tanglegram.pdf", width = 15, height = height) 

bac71_SCG_PAN_TANGLEGRAM %>% plot(common_subtrees_color_lines=FALSE,
                                  highlight_distinct_edges=FALSE,
                                  common_subtrees_color_branches=FALSE,
                                  highlight_branches_lwd=FALSE,
                                  lwd=3, 
                                  lab.cex = 0.65, edge.lwd = 2, 
                                  margin_inner = 20, 
                                  columns_width = c(1, 0.5, 1), 
                                  axes=FALSE,
                                  color_lines = line_colors3)



dev.off()


```

Print ANI tree

```{r, eval=FALSE}



# load ANI tree
df_newick <-  ape::read.tree("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/ANI/G_0213/ANIb_percentage_identity.newick")

  
# create dedrogram object
dend_df_newick <- ape::chronos(df_newick)
dend_df_newick_2 <- as.dendrogram.phylo(dend_df_newick)

# create rectangular lines object
ddata <- dendro_data(dend_df_newick_2, type = "rectangle")

ANI_dend_plot <- ggplot(segment(ddata)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), linewidth = 1) + 
  coord_flip() + 
  scale_y_reverse(expand = c(0.2, 0)) +
  scale_x_reverse() +
  theme_classic() +
  theme(axis.ticks.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y =  element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  theme(panel.spacing.x=unit(0.05, "lines"),
        plot.margin=unit(c(0,0,0,0), "cm"))

# save dend plot
ggsave("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/ANI/G_0213/ANIb_percentage_identity_dendrogram.pdf",ANI_dend_plot, width = 5, height = 7)



```



# 8. Grouping N. subflava and N. mucosa genomes based on ANI and Phylogenomics

We used a multifaceted approach to identify sub-species level genetic clades, hereafter referred to as subgroups, within the N. subflava and the N. mucosa major clades. Our methods included assessing Average Nucleotide Identity (ANI) and two single-copy core gene phylogenies: one based on the universal 71 bacterial single-copy core genes (SCCGs) and the other on SCCGs extracted from the pangenome for the respective focal genomes. For each species, we constructed three tanglegrams: one comparing the ANI-based hierarchical clustering tree with the 71 bacterial SCCGs, one comparing the ANI-based hierarchical clustering tree with the pangenome SCCGs, and one comparing the 71 bacterial SCCGs with the pangenome SCCGs. For ANI, we used the Gap statistic in conjunction with clusGap in the R package “cluster” v2.1.4 (34), which calculates a goodness-of-clustering measure (the “gap” statistic) for each possible number of k clusters. We also used the fviz_nbclust function in the R package “factoextra” v1.0.7 (35) to determine and visualize the optimal number of clusters using within-cluster sums of squares for each possible number of k clusters.

Please see the file *Grouping_N_sub_and_N_muc_genomes.Rmd* for relevant code used to identify genomic groups for N. subflava and N. mucosa genomes. 


# 9. Competative mapping metagenomic reads to selected set of reference genomes

##### 9.1. List of metagenome samples.

```{bash, eval=FALSE}

scp -r /Users/home/SPECIES_LEVEL_PANGENOMES/HMP_metadata/Final_1297_HMP_metagenomes_metadata.csv BARHAL_SERVER:/workspace/jmarkwelchlab/HMP_SAMPLE_DATA/Final_1297_HMP_metagenomes_metadata.csv

```

##### 9.2. Create symbolic links for HMP metagenomes. 

```{bash, eval=FALSE}

DIR_HMP_INFO=/workspace/jmarkwelchlab/HMP_SAMPLE_DATA
DIR_SITE=/storage3/data/02_QUALITY_FILTERED_HMP_SAMPLES
DIR_READS=/workspace/jmarkwelchlab/P_0003_Neisseriaceae/00_READS


for metagenome in `awk -F"," 'NR>1{print $3}' "$DIR_HMP_INFO/Final_1297_HMP_metagenomes_metadata.csv"`
do
    Internal_ID=$(awk -F"," -v meta="$metagenome" '$3 == meta {gsub(/"/, "", $1); print $1}' "$DIR_HMP_INFO/Final_1297_HMP_metagenomes_metadata.csv")
    
    metagenome_ID=$(awk -F"," -v meta="$metagenome" '$3 == meta {gsub(/"/, "", $3); print $3}' "$DIR_HMP_INFO/Final_1297_HMP_metagenomes_metadata.csv")
    
    ln -s "$DIR_SITE/${metagenome_ID}-QUALITY_PASSED_R1.fastq.gz" "$DIR_READS/${Internal_ID}_R1.fastq.gz"
    ln -s "$DIR_SITE/${metagenome_ID}-QUALITY_PASSED_R2.fastq.gz" "$DIR_READS/${Internal_ID}_R2.fastq.gz"
done


# check a link
ls -l /workspace/jmarkwelchlab/P_0003_Neisseriaceae/00_READS/PP_HC_HMP_S0054_01_R2.fastq.gz

```


Set up samples_id-QC_IDs.txt file, which is simply a header-less column containing sample IDs that match the prefixes of the _R1.fastq.gz and _R2.fastq.gz files. Make sure samples_id-QC_IDs.txt is in the main directory. 

```{bash, eval=FALSE}

projectID="P_0003_Neisseriaceae"
mainDIR="/workspace/jmarkwelchlab"
DIR_HMP_INFO=$mainDIR/HMP_SAMPLE_DATA

# copy samples_id-QC_IDs_686.txt to the samples_id-QC_IDs.txt file in the main directory for mapping
awk -F","  'NR>1 {gsub(/"/, "", $1); print $1}' "$DIR_HMP_INFO/Final_1297_HMP_metagenomes_metadata.csv" > $mainDIR/$projectID/samples_id-QC_IDs.txt

```

Directories that are needed to be cleared (i.e., move all stuff that is already in there from a previous into new directory named TEST_RUN05: 
05_MAPPING
06_SINGLE_PROFILE
07_MERGED_PROFILE
08_PROFILE_SUMMARY
13_DETECTED_GENOMES

copy scripts for MAPPING in parallel, profiling mapped reads, merging profiles my oral sites, and decomposing and summarizing merged profiles:

```{bash, eval=FALSE}
cp -r /workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/SCRIPTS  /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/Mapping_profiling_scripts

##### edit scripts
# Make sure min contig size 300
# Make sure internal script directories are ok (e.g., ...P_0003_Neisseriaceae/SCRIPTS/Mapping_profiling_scripts...)

# Single profiles
nano /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/Mapping_profiling_scripts/script-03-single_profiling_parallel.sh
/workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/Mapping_profiling_scripts/script-single_profile.sh

# Merging profiles
nano /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/Mapping_profiling_scripts/script-04-merged_profile_parallel.sh
/workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/Mapping_profiling_scripts/script-merged_profile-HP.sh
 /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/Mapping_profiling_scripts/script-merged_profile.sh
 
# decompose merged profile, estimate SCG taxonomy and summarize
nano /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/Mapping_profiling_scripts/script-05-decompose_summarize_profile_parallel.sh
/workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/Mapping_profiling_scripts/script-decompose_summarize_profile.sh

```


##### 9.3. Run mapping script.

```{bash, eval=FALSE}

projectID="P_0003_Neisseriaceae"

clusterize -n 15 -m jgiacomini@forsyth.org -log LOGS/mapping.log /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/Mapping_profiling_scripts/script-02-indexing_and_mapping_parallel.sh $projectID 

cat LOGS/mapping.log | grep "overall alignment rate" | wc -l

```


# 10. Profiling mapped reads 

##### 10.1 Run single profile script

```{bash, eval=FALSE}
projectID="P_0003_Neisseriaceae"

clusterize -n 15 -m jgiacomini@forsyth.org -log LOGS/profiling.log /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/Mapping_profiling_scripts/script-03-single_profiling_parallel.sh $projectID
```

Check results...
```{bash, eval=FALSE}
cat LOGS/profiling.log | grep "Happy" | wc -l
```

1297 samples and only 1283 occurrences of "Happy" for single profiling, indicating that 14 samples possibly failed. Also, running the merge profile script resulted in the following errors for some oral sites:

Config Error: The database '06_SINGLE_PROFILE/BM_HC_HMP_S0034_02/PROFILE.db' has no 'db_type'
              row in 'self'

Config Error: There is something wrong with your input databases. The group name 'default'
              should be common to all of them, but it doesn't seem to be the case :/ How did
              you end up with an anvi'o single profile database that doesn't have the
              'default' group in its additional layer data table? It is very likely that your
              profiling step failed for some reason for one or more of your databases :(
              

Unfortunately, Anvio did not tell me which single profile database caused the error/s. Fortunately there is a runlog asscocaiated with each sample in the profiling step. The following code will parse the information from each single profile run log. It will print the sample name and then whether or not the profile completed,  "happy", or did not "misisng". If the runlog is missing it will also print the sample name and the string "missing". 

```{bash, eval=FALSE}

for sample in `cat samples_id-QC_IDs.txt`
do
    # Print the sample name
    echo -n "$sample " >> Single_profile_runlog_Happy_test.txt
    
    # Check for the string "Happy" in the RUNLOG.txt
    if grep -q "Happy" 06_SINGLE_PROFILE/"$sample"/RUNLOG.txt; then
        # If "Happy" is found, print the line containing "Happy"
        echo "happy" >> Single_profile_runlog_Happy_test.txt
    else
        # If "Happy" is not found, print "missing"
        echo "missing" >> Single_profile_runlog_Happy_test.txt
    fi
done

cat Single_profile_runlog_Happy_test.txt | grep "missing" | wc -l 

```
        

Re-run single profiles for the 14 samples identified above:

```{bash, eval=FALSE}

# make list of missing samples
cat Single_profile_runlog_Happy_test.txt | grep "missing" | awk '{print $1}' > missing_samples_id-QC_IDs.txt

# move old single profiles to temp diretcory
mkdir FAILED_SINGLE_PROFILES
for sample in `cat  missing_samples_id-QC_IDs.txt`
do
mv 06_SINGLE_PROFILE/$sample FAILED_SINGLE_PROFILES/$sample
done

# Two samples missing entirely from 06_SINGLE_PROFILE directory:
# mv: cannot stat '06_SINGLE_PROFILE/PP_HC_HMP_S0135_02': No such file or directory
# mv: cannot stat '06_SINGLE_PROFILE/TD_HC_HMP_S0116_01': No such file or directory


# run the following script to profile the 14 missing samples
projectID="P_0003_Neisseriaceae"

clusterize -n 15 -m jgiacomini@forsyth.org -log LOGS/missing_samples_profiling.log /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/Mapping_profiling_scripts/script-03-single_profiling_parallel_missing_samples.sh $projectID
```

Check the new profiles:

```{bash}

for sample in `cat missing_samples_id-QC_IDs.txt`
do
    # Print the sample name
    echo -n "$sample " >> Missing_single_profile_runlog_Happy_test.txt
    
    # Check for the string "Happy" in the RUNLOG.txt
    if grep -q "Happy" 06_SINGLE_PROFILE/"$sample"/RUNLOG.txt; then
        # If "Happy" is found, print the line containing "Happy"
        echo "happy" >> Missing_single_profile_runlog_Happy_test.txt
    else
        # If "Happy" is not found, print "missing"
        echo "missing" >> Missing_single_profile_runlog_Happy_test.txt
    fi
done

cat Missing_single_profile_runlog_Happy_test.txt | grep "happy" | wc -l # 14
cat Missing_single_profile_runlog_Happy_test.txt | grep "missing" | wc -l # 0
```

All showed "happy"... success!!!
Also, we now see 1297 single profiles in the 06_SINGLE_PROFILE/ directory, yay! At least in terms of their run logs. That is not a guarantee that the data bases are ok though...foreshadowing...


##### 10.2 Merge profiles by oral site

Need to amend the merge profile for single site HP. I just manually copied the HP single profile directory into the merged profile directory.  

Note that I added the flags --overwrite-output-destinations & --skip-hierarchical-clustering 
```{bash, eval=FALSE}
projectID="P_0003_Neisseriaceae"

mkdir 07_MERGED_PROFILE/${projectID}_HP 
cp 06_SINGLE_PROFILE/HP_HC_HMP_S0195_01/* 07_MERGED_PROFILE/${projectID}_HP/

clusterize -n 15 -m jgiacomini@forsyth.org -log LOGS/merge_profiles.log /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/Mapping_profiling_scripts/script-04-merged_profile_parallel.sh $projectID  
```

Check the run logs for each merged profile. We want to see a "Happy" next to each oral site:

```{bash}
for site in TD BM HP KG PB PP PT SA TH; do

  status=$(cat 07_MERGED_PROFILE/P_0003_Neisseriaceae_${site}/RUNLOG.txt | grep -o "Happy")
  echo -e "$site\t$status"
  
done
```


##### 10.3 Summarize merged prfiles

The following script imports a collection to each oral site profile database that assigns a genome ID to the splits. It then summarizes the coverages for that collection, giving us coverages for each genome for each sample. 

```{bash, eval=FALSE}
projectID="P_0003_Neisseriaceae"

clusterize -n 20 -m jgiacomini@forsyth.org -log LOGS/decompose_summarize_2.log /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/Mapping_profiling_scripts/script-05-decompose_summarize_profile_parallel.sh $projectID  

```

Unfortunately I learned the hard way that anvi-summarize causes Barhal to crash when running the summarize_profile_parallel script. anvi-summarize uses a significant amount of memory. And so, we are still missing  summaries for two oral sites - TD and PP. See section "6.5 Check summaries" where I counted the number of rows and columns of the files that we need, which showed that all sites except for TD and PP were complete.  

In additon to running anvi-summarize, we also ran anvi-import-collection for each merged profile. The following command shows that both the TD and PP merged profile databases contain the required Genome collection and so we do not need to import the genome collection again.

```{bash, eval=FALSE}

anvi-show-collections-and-bins -p 07_MERGED_PROFILE/P_0003_Neisseriaceae_TD/PROFILE.db
anvi-show-collections-and-bins -p 07_MERGED_PROFILE/P_0003_Neisseriaceae_PP/PROFILE.db

```


The following scripts will run anvi-summarize for TD and PP, seperately. 

```{bash, eval=FALSE}

projectID="P_0003_Neisseriaceae"
clusterize -n 10 -m jgiacomini@forsyth.org -log LOGS/decompose_summarize_TD.log /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/Mapping_profiling_scripts/script-05-decompose_summarize_profile_NON_parallel.sh $projectID TD

projectID="P_0003_Neisseriaceae"
clusterize -n 20 -m jgiacomini@forsyth.org -log LOGS/decompose_summarize_PP.log /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/Mapping_profiling_scripts/script-05-decompose_summarize_profile_NON_parallel.sh $projectID PP


```


##### 10.4 Check summaries

```{bash, eval=FALSE}

declare -a siteArray=("BM" "TD" "PP" "PT" "TH" "KG" "PB" "SA" "HP")

for site in ${siteArray[@]}
do
    echo -e "$site"
    cat "08_PROFILE_SUMMARY/P_0003_Neisseriaceae_${site}-profile/bins_summary.txt" | wc -l
done



for site in ${siteArray[@]}
do
    echo -e "site\tdetection\tmean_coverage\tmean_coverage_Q2Q3"
    
   detection=$(awk -F'\t' '{print NF; exit}' "08_PROFILE_SUMMARY/P_0003_Neisseriaceae_${site}-profile/bins_across_samples/detection.txt")
   
    mean_coverage=$(awk -F'\t' '{print NF; exit}' "08_PROFILE_SUMMARY/P_0003_Neisseriaceae_${site}-profile/bins_across_samples/mean_coverage.txt")
    
    mean_coverage_Q2Q3=$(awk -F'\t' '{print NF; exit}' "08_PROFILE_SUMMARY/P_0003_Neisseriaceae_${site}-profile/bins_across_samples/mean_coverage_Q2Q3.txt")
    
    echo -e "$site\t$detection\t$mean_coverage\t$mean_coverage_Q2Q3"
done


# Number of rows for each should be 214
for site in ${siteArray[@]}
do
    echo -e "site\tdetection\tmean_coverage\tmean_coverage_Q2Q3"
    
   detection=$(cat "08_PROFILE_SUMMARY/P_0003_Neisseriaceae_${site}-profile/bins_across_samples/detection.txt" | wc -l)
   
    mean_coverage=$(cat "08_PROFILE_SUMMARY/P_0003_Neisseriaceae_${site}-profile/bins_across_samples/mean_coverage.txt" | wc -l)
    
    mean_coverage_Q2Q3=$(cat "08_PROFILE_SUMMARY/P_0003_Neisseriaceae_${site}-profile/bins_across_samples/mean_coverage_Q2Q3.txt" | wc -l)
    
    echo -e "$site\t$detection\t$mean_coverage\t$mean_coverage_Q2Q3"
done

```


# 11. Extract Mean depth and breadth of coverage data

Concatenate mean coverage (depth) and detection (breadth) for all oral sites
```{bash, eval=FALSE}

# set up variables
projectID="P_0003_Neisseriaceae"
mainDIR="/workspace/jmarkwelchlab"
DIR_DetectionGENOMES=$mainDIR/$projectID/13_DETECTED_GENOMES
DIR_SummaryPROF=$mainDIR/$projectID/08_PROFILE_SUMMARY

# Run for loop to extract and transpose data frames
for data in detection mean_coverage mean_coverage_Q2Q3
do
#variables
ref_site=SA
detTableTrans=$DIR_DetectionGENOMES/${projectID}-${data}.transposed.txt
detTableTransHeadless=$DIR_DetectionGENOMES/${projectID}-${data}.transposed.noheader.txt
detTSV=$DIR_DetectionGENOMES/${projectID}-${data}.txt
# transpose detection-tables to combine values from 9 oral sites
for site in PP PB BM KG TD PT TH HP SA 
do
projSite=${projectID}_${site}
siteDetTSV=$DIR_SummaryPROF/${projSite}-profile/bins_across_samples/${data}.txt
siteDetTSVtrans=$DIR_DetectionGENOMES/${projSite}-${data}.transposed.txt
anvi-script-transpose-matrix -o $siteDetTSVtrans $siteDetTSV
awk 'BEGIN{FS=OFS="\t"}NR>1{print $0}' $siteDetTSVtrans >> $detTableTransHeadless
done
# add header to detection-table (I choose ${ref_site}), transpose table, remove intermediate files
awk 'BEGIN{FS=OFS="\t"}NR==1{print $0}' $DIR_DetectionGENOMES/${projectID}_${ref_site}-${data}.transposed.txt > $detTableTrans
cat $detTableTransHeadless >> $detTableTrans
rm $DIR_DetectionGENOMES/${projectID}_*-${data}*txt $DIR_DetectionGENOMES/${projectID}-${data}*noheader*
anvi-script-transpose-matrix -o $detTSV $detTableTrans
rm $detTableTrans
done

```


Check the data - count the number of columns

P_0003_Neisseriaceae-detection.txt  P_0003_Neisseriaceae-mean_coverage_Q2Q3.txt  P_0003_Neisseriaceae-mean_coverage.txt

```{bash, eval=FALSE}


detection=$(awk -F'\t' '{print NF; exit}' "13_DETECTED_GENOMES/P_0003_Neisseriaceae-detection.txt")
mean_coverage=$(awk -F'\t' '{print NF; exit}' "13_DETECTED_GENOMES/P_0003_Neisseriaceae-mean_coverage.txt")
mean_coverage_Q2Q3=$(awk -F'\t' '{print NF; exit}' "13_DETECTED_GENOMES/P_0003_Neisseriaceae-mean_coverage_Q2Q3.txt")
    
echo -e "detection\tmean_coverage\tmean_coverage_Q2Q3"
echo -e "$detection\t$mean_coverage\t$mean_coverage_Q2Q3"

```

##### 11.1 Export total reads mapped data

Anvio reports mapping stats for  each sample, including the number of SNVs, SCVs, and indels. We can loop through each oral site profile and extract the info. However, we need to fix an issue where each oral site data frame has a different header order. The following R script will loop through the input files for the oral sites (BM TD PP PT PB TH KG SA HP), reorder columns so that they match. 

```{r, eval=FALSE}

# First need to arrange each data frame into the same order of columns
library("dplyr")

# Define the desired column order
new_cols <- c("layers", "total_reads_mapped", "total_reads_kept", "num_SNVs_reported", "num_SCVs_reported", "num_INDELs_reported")

# Loop through the input files and reorder columns BM TD PP PT PB TH KG SA HP 
for (file_name in c("/workspace/jmarkwelchlab/P_0003_Neisseriaceae/08_PROFILE_SUMMARY/P_0003_Neisseriaceae_TD-profile/misc_data_layers/default.txt",
                    "/workspace/jmarkwelchlab/P_0003_Neisseriaceae/08_PROFILE_SUMMARY/P_0003_Neisseriaceae_BM-profile/misc_data_layers/default.txt",
                    "/workspace/jmarkwelchlab/P_0003_Neisseriaceae/08_PROFILE_SUMMARY/P_0003_Neisseriaceae_PP-profile/misc_data_layers/default.txt",
                    "/workspace/jmarkwelchlab/P_0003_Neisseriaceae/08_PROFILE_SUMMARY/P_0003_Neisseriaceae_PT-profile/misc_data_layers/default.txt",
                    "/workspace/jmarkwelchlab/P_0003_Neisseriaceae/08_PROFILE_SUMMARY/P_0003_Neisseriaceae_PB-profile/misc_data_layers/default.txt",
                    "/workspace/jmarkwelchlab/P_0003_Neisseriaceae/08_PROFILE_SUMMARY/P_0003_Neisseriaceae_TH-profile/misc_data_layers/default.txt",
                    "/workspace/jmarkwelchlab/P_0003_Neisseriaceae/08_PROFILE_SUMMARY/P_0003_Neisseriaceae_KG-profile/misc_data_layers/default.txt",
                    "/workspace/jmarkwelchlab/P_0003_Neisseriaceae/08_PROFILE_SUMMARY/P_0003_Neisseriaceae_SA-profile/misc_data_layers/default.txt",
                    "/workspace/jmarkwelchlab/P_0003_Neisseriaceae/08_PROFILE_SUMMARY/P_0003_Neisseriaceae_HP-profile/misc_data_layers/default.txt"
                    )) {
  # Read in the data
  df <- read.delim(file_name, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  
  # Reorder the columns based on the desired order
  df <- df[, new_cols]
  
  # Write out the updated data to a new file
  write.table(df, file = paste0(file_name, "_new.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
}

```


Now we can combine files that we re-ordered above.
```{bash, eval=FALSE}

# set up variables
projectID="P_0003_Neisseriaceae"
mainDIR="/workspace/jmarkwelchlab"
DIR_DetectionGENOMES=$mainDIR/$projectID/13_DETECTED_GENOMES
DIR_SummaryPROF=$mainDIR/$projectID/08_PROFILE_SUMMARY
MAPPING_STATS=$DIR_DetectionGENOMES/${projectID}-mapping_stats_new.txt

# Run for loop to extract and transpose data frames

# add header
echo -e "layers\ttotal_reads_mapped\ttotal_reads_kept\tnum_SNVs_reported\tnum_SCVs_reported\tnum_INDELs_reported" >> $MAPPING_STATS

for site in BM TD PP PT PB TH KG SA HP 
do
MAPPING_STATS_SITE=$DIR_SummaryPROF/${projectID}_${site}-profile/misc_data_layers/default.txt_new.txt 
# add data
awk 'BEGIN{FS=OFS="\t"}NR>1{print $0}' $MAPPING_STATS_SITE >> $MAPPING_STATS
done
```

Download data from remote server.
```{bash, eval=FALSE}

scp -r BARHAL_SERVER:/workspace/jmarkwelchlab/P_0003_Neisseriaceae/13_DETECTED_GENOMES /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/13_DETECTED_GENOMES

```

##### 11.2 Sum coverage per species group

```{r}

library(tidyr)
library("dplyr")

Q2Q3_coverage_original <- read.csv("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/13_DETECTED_GENOMES/Strain_level_updated_Q2Q3_mean_coverage_data.csv", header = TRUE)

# select relevant columns
Q2Q3_coverage_original_selected <- Q2Q3_coverage_original %>% 
  dplyr::select(layers, bins, value)

# pivot wider
Q2Q3_coverage_original_selected_wide <- Q2Q3_coverage_original_selected %>%  
  pivot_wider(names_from = bins, values_from = value)

# mutate variable so that there is only site ID

Q2Q3_coverage_original_selected_wide_mutated <- Q2Q3_coverage_original_selected_wide %>% 
  mutate(site = case_when(grepl('TD_', layers) ~ 'TD', 
                                      grepl('KG_', layers) ~ 'KG', 
                                      grepl('PP_', layers) ~ 'SUPP',
                                      grepl('PB_', layers) ~ 'SUBP',
                                      grepl('SA_', layers) ~ 'SV',
                                      grepl('TH_', layers) ~ 'TH',
                                      grepl('PT_', layers) ~ 'PT',
                                      grepl('HP_', layers) ~ 'HP',
                                      grepl('BM_', layers) ~ 'BM'))

Q2Q3_coverage_original_selected_wide_mutated <- Q2Q3_coverage_original_selected_wide_mutated %>% dplyr::select(site, 2:(ncol(Q2Q3_coverage_original_selected_wide_mutated) - 1))
  
# Create file that contains list of genome IDs in first column and mean relative abundances for each site in other columns???
# Pivot the data to longer format
long_df <- pivot_longer(Q2Q3_coverage_original_selected_wide_mutated, cols = -site, names_to = "genomeIDs", values_to = "value")

long_df_with_groups <- merge(long_df, Final_Groups, by.x = "genomeIDs", by.y = "Pangenome_ID")

# Calculate the average relative abundance for each Genome_ID at each site
avg_Q2Q3_coverage_df_per_group <- long_df_with_groups %>%
  group_by(site, Group_ID) %>%
  summarise(avg_Q2Q3_coverage = mean(value, na.rm = TRUE)) %>%
  pivot_wider(names_from = site, values_from = avg_Q2Q3_coverage)

avg_Q2Q3_coverage_df_per_group <- avg_Q2Q3_coverage_df_per_group %>% 
  rename(groupIDs = Group_ID) 

# Save file
write.table(avg_Q2Q3_coverage_df_per_group, "/Users/home/SPECIES_LEVEL_PANGENOMES/TESTS/Classifying_genomes_to_oral_sites/Groups/With_HP/Neisseriaceae_mean_Q2Q3_coverage_per_group.txt", quote = FALSE, row.names = FALSE, sep = "\t")


# Remove HP site
avg_Q2Q3_coverage_df_per_group_no_HP <- avg_Q2Q3_coverage_df_per_group %>% 
  dplyr::select(-HP)


# Save file
write.table(avg_Q2Q3_coverage_df_per_group_no_HP, "/Users/home/SPECIES_LEVEL_PANGENOMES/TESTS/Classifying_genomes_to_oral_sites/Groups/Without_HP/Neisseriaceae_mean_Q2Q3_coverage_per_group_no_HP.txt", quote = FALSE, row.names = FALSE, sep = "\t")


########## Now sum coverages per group ########## 
Q2Q3_coverage_original_with_groups <- merge(Q2Q3_coverage_original_selected, Final_Groups, by.x = "bins", by.y = "Pangenome_ID")

sum_Q2Q3_coverage_df_per_group <- Q2Q3_coverage_original_with_groups %>%
  group_by(layers, Group_ID) %>%
  summarise(sum_Q2Q3_coverage = sum(value, na.rm = TRUE)) 

sum_Q2Q3_coverage_df_per_group <- sum_Q2Q3_coverage_df_per_group %>% 
  mutate(site = case_when(grepl('TD_', layers) ~ 'TD', 
                                      grepl('KG_', layers) ~ 'KG', 
                                      grepl('PP_', layers) ~ 'SUPP',
                                      grepl('PB_', layers) ~ 'SUBP',
                                      grepl('SA_', layers) ~ 'SV',
                                      grepl('TH_', layers) ~ 'TH',
                                      grepl('PT_', layers) ~ 'PT',
                                      grepl('HP_', layers) ~ 'HP',
                                      grepl('BM_', layers) ~ 'BM'))

# pivot wider
sum_Q2Q3_coverage_df_per_group_wide <- sum_Q2Q3_coverage_df_per_group %>%  
  dplyr::select(site, Group_ID, sum_Q2Q3_coverage) %>% 
  droplevels() %>% 
  pivot_wider(names_from = Group_ID, values_from = sum_Q2Q3_coverage) 

sum_Q2Q3_coverage_df_per_group_wide <- sum_Q2Q3_coverage_df_per_group_wide[,-1]


write.table(sum_Q2Q3_coverage_df_per_group_wide, "/Users/home/SPECIES_LEVEL_PANGENOMES/TESTS/Classifying_genomes_to_oral_sites/Groups/With_HP/Neisseriaceae_sum_Q2Q3_coverage_per_group.txt", quote = FALSE, row.names = FALSE, sep = "\t")
  
# Drop all rows that contain "HP" in "site: column

sum_Q2Q3_coverage_df_per_group_wide_no_HP <- sum_Q2Q3_coverage_df_per_group_wide %>% 
  filter(site != "HP")

write.table(sum_Q2Q3_coverage_df_per_group_wide_no_HP, "/Users/home/SPECIES_LEVEL_PANGENOMES/TESTS/Classifying_genomes_to_oral_sites/Groups/Without_HP/Neisseriaceae_sum_Q2Q3_coverage_per_group_no_HP.txt", quote = FALSE, row.names = FALSE, sep = "\t")

```

##### 11.3 Accumulated detection per species group

```{r}

library(tidyr)
library("dplyr")

Breadth_coverage_original <- read.csv("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/13_DETECTED_GENOMES/Strain_level_updated_detection_data.csv", header = TRUE)

Final_Groups <- read.table("/Users/home/SPECIES_LEVEL_PANGENOMES/TESTS/Classifying_genomes_to_oral_sites/Neisseriaceae_genome_groups.txt", header = TRUE, sep = "\t")
 
# select relevant columns
Detection_original_selected <- Breadth_coverage_original %>% 
  dplyr::select(layers, bins, detection_binary)

# Merge with species groups info
Detection_original_with_groups <- merge(Detection_original_selected, Final_Groups, by.x = "bins", by.y = "Pangenome_ID")

# Sum detection and mutate detection so that values >=1 become binary
sum_detection_df_per_group <- Detection_original_with_groups %>%
  group_by(layers, Group_ID) %>%
  summarise(sum_detection = sum(detection_binary)) %>% 
  mutate(accumulated_detection = ifelse(sum_detection >= 1, 1, 0))


# Add site column
sum_detection_df_per_group <- sum_detection_df_per_group %>% 
  mutate(site = case_when(grepl('TD_', layers) ~ 'TD', 
                                      grepl('KG_', layers) ~ 'KG', 
                                      grepl('PP_', layers) ~ 'SUPP',
                                      grepl('PB_', layers) ~ 'SUBP',
                                      grepl('SA_', layers) ~ 'SV',
                                      grepl('TH_', layers) ~ 'TH',
                                      grepl('PT_', layers) ~ 'PT',
                                      grepl('HP_', layers) ~ 'HP',
                                      grepl('BM_', layers) ~ 'BM'))

# pivot wider 
sum_detection_df_per_group_wide <- sum_detection_df_per_group %>%  
  dplyr::select(site, Group_ID, accumulated_detection) %>% 
  droplevels() %>% 
  pivot_wider(names_from = Group_ID, values_from = accumulated_detection) 

sum_detection_df_per_group_wide <- sum_detection_df_per_group_wide[,-1]

# Save wide df

write.table(sum_detection_df_per_group_wide, "/Users/home/SPECIES_LEVEL_PANGENOMES/TESTS/Classifying_genomes_to_oral_sites/Groups/With_HP/Neisseriaceae_accum_detection_per_group.txt", quote = FALSE, row.names = FALSE, sep = "\t")
  
# Drop all rows that contain "HP" in "site: column
sum_detection_df_per_group_wide_no_HP <- sum_detection_df_per_group_wide %>% 
  filter(site != "HP")

write.table(sum_detection_df_per_group_wide_no_HP, "/Users/home/SPECIES_LEVEL_PANGENOMES/TESTS/Classifying_genomes_to_oral_sites/Groups/Without_HP/Neisseriaceae_accum_detection_per_group_no_HP.txt", quote = FALSE, row.names = FALSE, sep = "\t")


# Pivot from wide to long form just incase needed
sum_detection_df_per_group_wide_no_HP_long <- pivot_longer(sum_detection_df_per_group_wide_no_HP, cols = -site, names_to = "groupIDs", values_to = "Detection_Binary") %>%
  rename(Site_ID = site) 

write.table(sum_detection_df_per_group_wide_no_HP_long, "/Users/home/SPECIES_LEVEL_PANGENOMES/TESTS/Classifying_genomes_to_oral_sites/Groups/Without_HP/Neisseriaceae_accum_detection_per_group_no_HP_Long.txt", quote = FALSE, row.names = FALSE, sep = "\t")

```

# 12. FIGURES: Detection, relative abundance and mean coverage

See *Neisseriaceae_detection_mean_coverage_rel_abundance_plots.Rmd* for detection, relative abundance and mean coverage plots.


# 13. Multi-site-preference classification

See *Neisseriaceae_Multi_site_preference_classification.Rmd*


# 14. Metapangenome

## 14.1. Generate oral site layers for metapangenome. 

Run after profile summary complete.

```{bash Generate oral site layers for metapangenome}
# LAYERS from metadata and merge with mapping reads and relative mapping
projectID="P_0003_Neisseriaceae"
mainDIR="/workspace/jmarkwelchlab"
DIR_Data="$mainDIR/$projectID/DATA"
DIR_SummaryPROF="$mainDIR/$projectID/08_PROFILE_SUMMARY"
mapIntermediate="$DIR_Data/mapping-reads.intermediate.tsv"
samplesMetadata="$DIR_Data/sequencing-metadata-count-QC.txt"
LAYER_ORDERS="$DIR_Data/$projectID-1297-layer_orders.txt"
sequencing_metadata="$mainDIR/HMP_SAMPLE_DATA/Final_1297_HMP_metagenomes_metadata.txt"
samples="$mainDIR/$projectID/samples_id-QC_IDs.txt"

# Update samples metadata
input="HMP_SAMPLE_DATA/Final_1297_HMP_metagenomes_metadata.csv"
output="HMP_SAMPLE_DATA/Final_1297_HMP_metagenomes_metadata.txt"
# Remove quotes and convert to tab-separated
sed 's/"//g' $input | awk 'BEGIN {FS=","; OFS="\t"} {$1=$1; print}' > $output

# Get sample IDs and total number of mapped reads
for site in PP PB BM KG TD PT HP TH SA 
do
siteLayers=$DIR_SummaryPROF/${projectID}_${site}-profile/misc_data_layers/default.txt
awk -F'\t' -v OFS="\t" 'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}} { print $(f["layers"]), $(f["total_reads_mapped"])}' $siteLayers | awk 'BEGIN{FS=OFS="\t"}NR>1{print $1,$2}' >> $mapIntermediate
done


# Write the header
echo -e "Sample_ID\tOral_site\tSubject_gender\tTotal_reads\tMapping_reads\tRel_map_reads" > $samplesMetadata
for metagenome in `cat samples_id-QC_IDs.txt`
do
META=$(cat $sequencing_metadata | grep -w "$metagenome" | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$13,$19}')
MAP=$(cat $mapIntermediate | grep -w "$metagenome" | awk 'BEGIN{FS=OFS="\t"}{print $2}')
echo -e "$META\t$MAP" | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6=((100 * $5 )/ ($4 * 2))}' >> $samplesMetadata
done
rm $mapIntermediate



# generate LAYER_ORDERS file
# order samples using multiple conditions (intermediate files)
for site in PP PB BM KG TD PT HP TH SA 
do
# Abundance
awk 'BEGIN{FS=OFS="\t"}NR>1{ print $1,$2,$6}' $samplesMetadata | sort -k2,2 -k3,3n | cut -f1 > Abundance.intermediate.layer_orders.tsv
# Site_abundance
awk -v SITE="$site" 'BEGIN{FS=OFS="\t"}NR>1{if( $1 ~SITE) print $1,$2,$6}' $samplesMetadata | sort -k2,2 -k3,3n | cut -f1 >> Site_abundance.intermediate.layer_orders.tsv
# Site_abundance_reverse
awk -v SITE="$site" 'BEGIN{FS=OFS="\t"}NR>1{if( $1 ~SITE) print $1,$2,$6}' $samplesMetadata | sort -k2,2 -k3,3nr | cut -f1 >> Site_abundance_reverse.intermediate.layer_orders.tsv
# Site_reverse_abundance
awk -v SITE="$site" 'BEGIN{FS=OFS="\t"}NR>1{if( $1 ~SITE) print $1,$2,$6}' $samplesMetadata | sort -k2,2r -k3,3n | cut -f1 >> Site_reverse_abundance.intermediate.layer_orders.tsv
# Site_reverse_abundance_reverse
awk -v SITE="$site" 'BEGIN{FS=OFS="\t"}NR>1{if( $1 ~SITE) print $1,$2,$6}' $samplesMetadata | sort -k2,2r -k3,3nr | cut -f1 >> Site_reverse_abundance_reverse.intermediate.layer_orders.tsv
# Site_mapped_reads
awk -v SITE="$site" 'BEGIN{FS=OFS="\t"}NR>1{if( $1 ~SITE) print $1,$2,$5}' $samplesMetadata | sort -k2,2 -k3,3n | cut -f1 >> Site_mapped_reads.intermediate.layer_orders.tsv
# Site_mapped_reads_reverse
awk -v SITE="$site" 'BEGIN{FS=OFS="\t"}NR>1{if( $1 ~SITE) print $1,$2,$5}' $samplesMetadata | sort -k2,2 -k3,3nr | cut -f1 >> Site_mapped_reads_reverse.intermediate.layer_orders.tsv
# Site_reverse_mapped_reads
awk -v SITE="$site" 'BEGIN{FS=OFS="\t"}NR>1{if( $1 ~SITE) print $1,$2,$5}' $samplesMetadata | sort -k2,2r -k3,3n | cut -f1 >> Site_reverse_mapped_reads.intermediate.layer_orders.tsv
# Site_reverse_mapped_reads_reverse
awk -v SITE="$site" 'BEGIN{FS=OFS="\t"}NR>1{if( $1 ~SITE) print $1,$2,$5}' $samplesMetadata | sort -k2,2r -k3,3nr | cut -f1 >> Site_reverse_mapped_reads_reverse.intermediate.layer_orders.tsv
# Site_total_reads
awk -v SITE="$site" 'BEGIN{FS=OFS="\t"}NR>1{if( $1 ~SITE) print $1,$2,$4}' $samplesMetadata | sort -k2,2 -k3,3n | cut -f1 >> Site_total_reads.intermediate.layer_orders.tsv
# Site_total_reads_reverse
awk -v SITE="$site" 'BEGIN{FS=OFS="\t"}NR>1{if( $1 ~SITE) print $1,$2,$4}' $samplesMetadata | sort -k2,2 -k3,3nr | cut -f1 >> Site_total_reads_reverse.intermediate.layer_orders.tsv
# Site_reverse_total_reads
awk -v SITE="$site" 'BEGIN{FS=OFS="\t"}NR>1{if( $1 ~SITE) print $1,$2,$4}' $samplesMetadata | sort -k2,2r -k3,3n | cut -f1 >> Site_reverse_total_reads.intermediate.layer_orders.tsv
# Site_reverse_total_reads_reverse
awk -v SITE="$site" 'BEGIN{FS=OFS="\t"}NR>1{if( $1 ~SITE) print $1,$2,$4}' $samplesMetadata | sort -k2,2r -k3,3nr | cut -f1 >> Site_reverse_total_reads_reverse.intermediate.layer_orders.tsv
done


# create LAYER ORDERS file
echo -e 'item_name\tdata_type\tdata_value' > $LAYER_ORDERS
# linearize order and append
for file in *intermediate.layer_orders.tsv
do
FILE_NAME=$(echo $file | awk -F'.' '{print $1}')
paste -sd"," $file | awk -v file_name="$FILE_NAME" 'BEGIN{FS=OFS="\t"}{print file_name,"basic",$0}' >> $LAYER_ORDERS
done


# transpose, linearize order and append
for file in *intermediate.layer_orders.tsv
do
FILE_NAME=$(echo $file | awk -F'.' '{print "Inverted_"$1}')
cat $file | tac | paste -sd","| awk -v file_name="$FILE_NAME" 'BEGIN{FS=OFS="\t"}{print file_name,"basic",$0}' >> $LAYER_ORDERS
done

# remove intermediate files
rm *intermediate.layer_orders.tsv

```


## 14.2. Build metpangenome.

```{bash, eval=FALSE}

projectID="P_0003_Neisseriaceae"

clusterize -n 10 -m jgiacomini@forsyth.org -log LOGS/metapangenome_summarise.log /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/script-06-metapangenome_summarize_pan_parallel.sh $projectID default

```

Need to repeat for PP and HP. 
```{bash}

#First remove HMP_PP and HMP_HP
rm -r 09_PANGENOME/P_0003_Neisseriaceae_G_0213-pangenome/HMP_PP 09_PANGENOME/P_0003_Neisseriaceae_G_0213-pangenome/HMP_HP

# Run again
clusterize -n 10 -m jgiacomini@forsyth.org -log LOGS/metapangenome_summarise_run_2.log /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/script-06-metapangenome_summarize_pan_parallel.sh $projectID default

```

## 14.3. Environmental core & accessory data.

```{bash ECG and EAG}

projectID="P_0003_Neisseriaceae"
mainDIR="/workspace/jmarkwelchlab"
DIR_Pangenome=$mainDIR/$projectID/09_PANGENOME/P_0003_Neisseriaceae_G_0213-pangenome
genomeDB=$DIR_Pangenome/${projectID}-GENOMES.db
panDB=$DIR_Pangenome/${projectID}-RESULTS/${projectID}-PAN.db
panLayersOrder=$DIR_Pangenome/${projectID}-layer_orders.txt

# combining metapangenome
# get environmental core/accessory data
for site in BM TD PP PT PB TH KG SA HP 
do
siteDir=$DIR_Pangenome/HMP_${site}
panDB=$siteDir/${projectID}-RESULTS/${projectID}-PAN.db
itemsTable=$siteDir/${projectID}-${site}-items.txt
anvi-export-misc-data -p $panDB -t items -o $itemsTable
cat $itemsTable | awk -F'\t' -v OFS="\t" 'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}} { print $(f["ECGs_and_EAGs!EAG"]),$(f["ECGs_and_EAGs!ECG"]),$(f["ECGs_and_EAGs!NA"]),$(f["EAG_ECG_ratio"])}' | awk -v SITE="${site}_" -v OFS="\t" 'NR==1{print SITE$1,SITE$2,SITE$3,SITE$4}NR>1{print $1,$2,$3,$4}' > ${itemsTable}.intermediate
done


######## CHECK ############################################################
for site in BM TD PP PT PB TH KG SA HP 
do
    check_metapan(){
    cat 09_PANGENOME/P_0003_Neisseriaceae_G_0213-pangenome/HMP_${site}/P_0003_Neisseriaceae-${site}-items.txt.intermediate | grep -q "ECGs_and_EAGs"
    }
    if ! check_metapan; then
        echo "${site} ***Incomplete***"
    else
        echo "${site} Completed"
    fi
done

```

BM Completed
TD Completed
PP Completed
PT Completed
PB Completed
TH Completed
KG Completed
SA Completed
HP Completed

Import gene cluster ECG and EAG data to the metapangenome.

```{bash, eval=FALSE}

projectID="P_0003_Neisseriaceae"
mainDIR="/workspace/jmarkwelchlab"
DIR_Pangenome=$mainDIR/$projectID/09_PANGENOME/P_0003_Neisseriaceae_G_0213-pangenome
panDB=$DIR_Pangenome/${projectID}-RESULTS/${projectID}-PAN.db
environmentalGenes=$DIR_Pangenome/environmental_genes.tsv

cat $DIR_Pangenome/HMP_BM/${projectID}-BM-items.txt | cut -f1 > $environmentalGenes.intermediate
paste $environmentalGenes.intermediate $DIR_Pangenome/HMP_*/${projectID}-*-items.txt.intermediate > $environmentalGenes
rm $environmentalGenes.intermediate

# import to pangenomeDB
anvi-import-misc-data -p $panDB -t items $environmentalGenes

```

# 15. Functional enrichment analyses

Metabolic capabilities of genomes were predicted using the Anvi'o script anvi-estimate-metabolism (with parameters: --kegg-output-modes modules). The metabolic pathways considered by this script are those outlined in the KEGG MODULE Database and defined by KEGG Orthologs (KOs) (28). Each KO represents a specific gene function, and a KEGG module is a collection of KOs that work together to complete the steps of a metabolic pathway. We used the default completion threshold of 0.75, which scores a metabolic pathway as “complete” within a genome when at least 75% of the enzymes in the pathway are present in the genome. We then used the Anvi’o script anvi-compute-metabolic-enrichment to identify complete metabolic pathways differentially enriched in one set of genomes compared to another based on their habitat preferences. 

To identify COG functional annotations that are differentially enriched or depleted in one set of genomes compared to another, we used the Anvi’o script anvi-compute-functional-enrichment. The script associated each gene cluster with the most frequently annotated function and generated a frequency table of functions across genomes. An enrichment test was then conducted using a generalized linear model with a logit linkage function to obtain the enrichment score and an adjusted p-value (q-value). For both the KEGG and COG analyses, we considered modules or functions to be significantly enriched if they had a q-value of less than 0.01.


Overall, three primary habitat preferences emerge among the Neisseriaceae taxa. These include taxa that predominantly prefer dental plaque (SUPP and SUBP), the tongue dorsum (TD), or the keratinized gingiva (KG). The majority of Neisseriaceae taxa specialize in dental plaque, with 11 out of the 14 species detected in the oral cavity exhibiting this preference. Only one taxon, N. subflava, is specialized for the tongue dorsum habitat, while two species, N. cinerea and S. muelleri, specialize in the keratinized gingiva.

#### 15.1. Compute pangenome functional enrichment

```{bash}

projectID="P_0003_Neisseriaceae"
mainDIR="/workspace/jmarkwelchlab"
PAN="$mainDIR/$projectID/09_PANGENOME/P_0003_Neisseriaceae_G_0213-pangenome/P_0003_Neisseriaceae-RESULTS/P_0003_Neisseriaceae-PAN.db"
GENOMES="$mainDIR/$projectID/09_PANGENOME/P_0003_Neisseriaceae_G_0213-pangenome/P_0003_Neisseriaceae-GENOMES.db"
OUT_DIR="$mainDIR/$projectID/20_FUNCTIONAL_ENRICHMENT"

# Run loop for annotation sources...
for annotation_source in Pfam KEGG_Class KEGG_BRITE COG20_PATHWAY KEGG_Module COG20_CATEGORY KOfam COG20_FUNCTION
do
anvi-compute-functional-enrichment-in-pan -p $PAN -g $GENOMES -o $OUT_DIR/${annotation_source}_Pangenome_Enrichment_GRANULAR.txt --category-variable Habitat_Preference_Granular --annotation-source $annotation_source --functional-occurrence-table-output $OUT_DIR/${annotation_source}_Pangenome_FUNC_OCCURRENCE_GRANULAR.TXT
done

# ran separately for fucntional enrichment of Gene cluster identity...
anvi-compute-functional-enrichment-in-pan -p $PAN -g $GENOMES -o $OUT_DIR/GC_IDENTITY_Pangenome_Enrichment_GRANULAR.txt --category-variable Habitat_Preference_Granular --functional-occurrence-table-output $OUT_DIR/GC_IDENTITY_Pangenome_FUNC_OCCURRENCE_GRANULAR.TXT --annotation-source IDENTITY --include-gc-identity-as-function

```

#### 15.2. Compute metabolic pathway enrichment

##### Set up grouping file for anvi-compute-metabolic-enrichment

```{r, eval=FALSE}

library("dplyr")

df <- read.table("/workspace/jmarkwelchlab/P_0003_Neisseriaceae/12_FUNCTIONAL_ANNOTATION/Refined_habitat_prefs.txt", header = TRUE, sep = "\t")

metabolic_enrichment_groups <- df %>% 
  dplyr::select(item, Habitat_Preference_Granular) %>% 
  dplyr::rename(name = item,
                group = Habitat_Preference_Granular)

write.table(metabolic_enrichment_groups, "/workspace/jmarkwelchlab/P_0003_Neisseriaceae/12_FUNCTIONAL_ANNOTATION/habitat_prefs_site_granular_version_KEGG.txt", row.names = FALSE, quote = FALSE, sep = "\t")

```

##### Run anvi-compute-metabolic-enrichment


```{bash, eval = FALSE}

projectID="P_0003_Neisseriaceae"
mainDIR="/workspace/jmarkwelchlab"
DIR_Annotation=$mainDIR/$projectID/12_FUNCTIONAL_ANNOTATION
OUT_DIR="$mainDIR/$projectID/20_FUNCTIONAL_ENRICHMENT"
GRANULAR_pangenome_groupings_for_KEGG="$DIR_Annotation/habitat_prefs_site_granular_version_KEGG.txt"

# Run ernichment test. Module completion 75% default
anvi-compute-metabolic-enrichment -M $DIR_Annotation/KEGG_Metabolic_pathways/${projectID}-MetabolicPaths_modules.txt -G $GRANULAR_pangenome_groupings_for_KEGG --module-completion-threshold 0.75 -o $OUT_DIR/KEGG_Metabolic_Pathway_Enrichment_completion_75_GRANULAR.txt

# Run ernichment test. Module completion 100% 
anvi-compute-metabolic-enrichment -M $DIR_Annotation/KEGG_Metabolic_pathways/${projectID}-MetabolicPaths_modules.txt -G $GRANULAR_pangenome_groupings_for_KEGG --module-completion-threshold 1.00 -o $OUT_DIR/KEGG_Metabolic_Pathway_Enrichment_completion_100_GRANULAR.txt
```

#### 15.3. Closely related organisms predominantly inhabit different environments
 
Variants of N. mucosa predominantly preferred SUPP and SUBP habitats. In contrast, all but one variant of N. subflava showed a preference for TD habitats, with some variants also favoring PT, TH, and SV environments.
 
We can compare the N. mucosa genomes to the N. subflava genomes that prefer TD to give us some insight into how closely related organisms can predominantly inhabit different environments. 

 
##### Make grouping file for comparisons

```{r, eval = FALSE}

library("dplyr")

df <- read.table("/workspace/jmarkwelchlab/P_0003_Neisseriaceae/12_FUNCTIONAL_ANNOTATION/Refined_habitat_prefs.txt", header = TRUE, sep = "\t")


 mucosa_vs_subflava <- df %>% 
  dplyr::select(item, Group_ID, Habitat_Preference_Granular) %>% 
  dplyr::rename(name = item,
                group = Habitat_Preference_Granular) %>% 
  dplyr::filter(Group_ID == "N_mucosa_I" |
                Group_ID == "N_mucosa_VIII" |
                Group_ID == "N_mucosa_IV" |
                Group_ID == "N_mucosa_V" |
                Group_ID == "N_mucosa_VI" |
                Group_ID == "N_mucosa_VII" |
                Group_ID == "N_subflava_I" |
                Group_ID == "N_subflava_II" |
                Group_ID == "N_subflava_III" |
                Group_ID == "N_subflava_IV" |
                Group_ID == "N_subflava_VII" |
                Group_ID == "N_subflava_VIII") %>% 
  dplyr::select(name, group) %>% 
  dplyr::mutate(group = case_when(grepl('TONGUE', group) ~ 'N_subflava',
                                  grepl('PLAQUE', group) ~ 'N_mucosa'))

write.table(mucosa_vs_subflava, "/workspace/jmarkwelchlab/P_0003_Neisseriaceae/12_FUNCTIONAL_ANNOTATION/mucosa_vs_subflava_groups.txt", row.names = FALSE, quote = FALSE, sep = "\t")


```

##### COG20 Functions

```{bash, eval=FALSE}

projectID="P_0003_Neisseriaceae"
mainDIR="/workspace/jmarkwelchlab"
DIR_Annotation="$mainDIR/$projectID/12_FUNCTIONAL_ANNOTATION"
DIR_Pangenome="$mainDIR/$projectID/09_PANGENOME"
INTERNAL="$DIR_Pangenome/${projectID}_G_0213-pangenome/internal_${projectID}.txt"
GENOMES="$DIR_Pangenome/${projectID}_G_0213-pangenome/${projectID}-GENOMES.db"
GROUP_file="$DIR_Annotation/mucosa_vs_subflava_groups.txt"

#mkdir $DIR_Annotation/COG20_Function

anvi-display-functions -g $GENOMES \
         --groups-txt $GROUP_file \
         --annotation-source COG20_FUNCTION \
         --profile-db $DIR_Annotation/COG20_Function/mucosa_vs_subflava-COG20_FUNCTION-PROFILE.db \
         | tee LOGS/mucosa_vs_subflava-COG20_FUNCTION.log


mv /usr/local/tmp/tmpdbjmnfkb/FUNC_OCCURENCE_STATS.txt $DIR_Annotation/COG20_Function/mucosa_vs_subflava-COG20_FUNCTION-FUNC_OCCURENCE_STATS.txt

mv /usr/local/tmp/tmpdbjmnfkb/FUNC_ENRICHMENT_OUTPUT.txt $DIR_Annotation/COG20_Function/mucosa_vs_subflava-COG20_FUNCTION-FUNC_ENRICHMENT_OUTPUT.txt

```         

##### Pfams


```{bash, eval=FALSE}

projectID="P_0003_Neisseriaceae"
mainDIR="/workspace/jmarkwelchlab"
DIR_Annotation="$mainDIR/$projectID/12_FUNCTIONAL_ANNOTATION"
DIR_Pangenome="$mainDIR/$projectID/09_PANGENOME"
INTERNAL="$DIR_Pangenome/${projectID}_G_0213-pangenome/internal_${projectID}.txt"
GENOMES="$DIR_Pangenome/${projectID}_G_0213-pangenome/${projectID}-GENOMES.db"
GROUP_file="$DIR_Annotation/mucosa_vs_subflava_groups.txt"

#mkdir $DIR_Annotation/Pfam

anvi-display-functions -g $GENOMES \
         --groups-txt $GROUP_file \
         --annotation-source Pfam \
         --profile-db $DIR_Annotation/COG20_Function/mucosa_vs_subflava-Pfam-PROFILE.db \
         | tee LOGS/mucosa_vs_subflava-Pfam


mv /usr/local/tmp/tmp_vvt9v58/FUNC_OCCURENCE_STATS.txt $DIR_Annotation/Pfam/mucosa_vs_subflava-Pfam-FUNC_OCCURENCE_STATS.txt

mv /usr/local/tmp/tmp_vvt9v58/FUNC_ENRICHMENT_OUTPUT.txt $DIR_Annotation/Pfam/mucosa_vs_subflava-Pfam-FUNC_ENRICHMENT_OUTPUT.txt

```  

##### KEGG metabolic pathway enrichment


```{bash, eval = FALSE}

projectID="P_0003_Neisseriaceae"
mainDIR="/workspace/jmarkwelchlab"
DIR_Annotation=$mainDIR/$projectID/12_FUNCTIONAL_ANNOTATION
OUT_DIR="$mainDIR/$projectID/20_FUNCTIONAL_ENRICHMENT"
GRANULAR_pangenome_groupings_for_KEGG="$DIR_Annotation/mucosa_vs_subflava_groups.txt"

# Run ernichment test. Module completion 75% default
anvi-compute-metabolic-enrichment -M $DIR_Annotation/KEGG_Metabolic_pathways/${projectID}-MetabolicPaths_modules.txt -G $GRANULAR_pangenome_groupings_for_KEGG --module-completion-threshold 0.75 -o $OUT_DIR/KEGG_Metabolic_Pathway_Enrichment_completion_75_GRANULAR_mucosa_vs_subflava.txt

# Run ernichment test. Module completion 100% 
anvi-compute-metabolic-enrichment -M $DIR_Annotation/KEGG_Metabolic_pathways/${projectID}-MetabolicPaths_modules.txt -G $GRANULAR_pangenome_groupings_for_KEGG --module-completion-threshold 1.00 -o $OUT_DIR/KEGG_Metabolic_Pathway_Enrichment_completion_100_GRANULAR_mucosa_vs_subflava.txt


```
 

#### 15.4. Send ALL FUNCTIONAL_ENRICHMENT results to local (RUN ON LOCAL)

```{bash, eval=FALSE}

scp -r BARHAL_SERVER:/workspace/jmarkwelchlab/P_0003_Neisseriaceae/20_FUNCTIONAL_ENRICHMENT /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/

```


# 16. Nitrogen metabolism

## 16.1 Occurence of denitrification genes 

Build a table that summarizes the occurence of nitrate reduction genes in the Neisseriacea pangenome. 

1. Get list of all genomes.
2. Get list of focal genes - COG20 acc. and KOfams acc.
3. Create R function to search through COG20 functions and KOfams files for occurrence of focal genes in each of the genome. Output 1 for yes occurrence and 0 for no occurrence. 

```{r, eval=FALSE}

COG20_occurrence <- read.table("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/20_FUNCTIONAL_ENRICHMENT/Neisseriaceae_pangenome_G0213_COG20_Functions.txt", header = TRUE, sep = "\t", fill = TRUE, quote = "")
  
KOfam_occurrence <- read.table("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/20_FUNCTIONAL_ENRICHMENT/Neisseriaceae_pangenome_G0213_KOfams.txt", header = TRUE, sep = "\t", fill = TRUE, quote = "")

Genomes <- read.table("/Users/home/SPECIES_LEVEL_PANGENOMES/TESTS/Classifying_genomes_to_oral_sites/Neisseriaceae_genome_groups.txt", header = TRUE, sep = "\t")

Gene_Info <- read.table("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/Nitrate/Nitrate_metabolism_gene_info.txt", header = TRUE, sep = "\t")


# Load necessary libraries
library("dplyr")
library(tidyr)

# Rename the function column in COG20_occurrence and KOfam_occurrence
COG20_occurrence <- COG20_occurrence %>% rename(gene_function = function.)
KOfam_occurrence <- KOfam_occurrence %>% rename(gene_function = function.)

# Define the function
generate_occurrence_table <- function(COG20_occurrence, KOfam_occurrence, Genomes, Gene_Info) {
  
  # Create a list of unique Pangenome_IDs
  pangenome_ids <- unique(Genomes$Pangenome_ID)
  
  # Initialize an empty data frame for the occurrence table
  occurrence_table <- data.frame()
  
  # Iterate over each Pangenome_ID
  for (pangenome_id in pangenome_ids) {
    
    # Subset the Genomes data for the current Pangenome_ID
    genome_subset <- Genomes %>% filter(Pangenome_ID == pangenome_id)
    
    # Get the corresponding Group_ID
    group_id <- unique(genome_subset$Group_ID)
    
    # Check for the presence of KOfam_acc in KOfam_occurrence
    kofam_present <- KOfam_occurrence %>%
      filter(genome_name %in% genome_subset$Pangenome_ID) %>%
      pull(accession)
    
    # Check for the presence of COG_acc in COG20_occurrence
    cog_present <- COG20_occurrence %>%
      filter(genome_name %in% genome_subset$Pangenome_ID) %>%
      pull(accession)
    
    # Initialize a data frame for the current Pangenome_ID with all zeros
    current_table <- data.frame(
      Pangenome_ID = pangenome_id,
      Group_ID = group_id,
      matrix(0, nrow = 1, ncol = nrow(Gene_Info))
    )
    colnames(current_table)[3:ncol(current_table)] <- Gene_Info$Gene_name
    
    # Update the occurrence table with the presence (1) of genes
    for (i in 1:nrow(Gene_Info)) {
      gene_name <- Gene_Info$Gene_name[i]
      kofam_acc <- Gene_Info$KOfam_acc[i]
      cog_acc <- Gene_Info$COG_acc[i]
      
      if (kofam_acc %in% kofam_present || cog_acc %in% cog_present) {
        current_table[1, gene_name] <- 1
      }
    }
    
    # Bind the current table to the occurrence table
    occurrence_table <- bind_rows(occurrence_table, current_table)
  }
  
  return(occurrence_table)
}

# Run function 
occurrence_table <- generate_occurrence_table(COG20_occurrence, KOfam_occurrence, Genomes, Gene_Info)

# Merge results with genome habitat preference results
strain_level_habitat_classification_res <- read.table("/Users/home/SPECIES_LEVEL_PANGENOMES/TESTS/Classifying_genomes_to_oral_sites/Groups/Without_HP/Final_results_strain_level.txt", header = TRUE, sep = "\t")

merged_res <- merge(occurrence_table, strain_level_habitat_classification_res, by = "Pangenome_ID")

# clean up

merged_res <- merged_res %>% 
  dplyr::select(-c(Group_ID.y)) %>% 
  dplyr::rename(Group_ID = Group_ID.x) %>% 
  dplyr::select(Pangenome_ID, Group_ID, ABUNDANT_IN, PREVALENT_IN, COMBINED_CLASSIFICATION, GENERAL_SITE_CLASS, narG,
                narH, narI, narJ, nirK, nirS, norB, norC, 
                nosZ, nosR, nosD, nosF, nosY, nosL, nosX,
                narQ, narL, narX, narP,fnr, nsrR, 
                nrt2, 
                napA, napB,nirB, nirD, nrfA, nrfH )

# Save the occurrence table
write.table(merged_res, "/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/Nitrate/occurrence_table.txt", sep = "\t", row.names = FALSE, quote = FALSE)



```


## 16.2. Coverage of denitrification genes 

1. Parse list of gene caller IDs for each denitrification gene.
2. Get coverages of selected gene caller IDs and make bar plots.

```{r, eval = FALSE}

# Define the function
generate_gene_callers_list <- function(COG20_occurrence, KOfam_occurrence, Genomes, Gene_Info) {
  
  # Merge COG20_occurrence and KOfam_occurrence
  occurrences <- bind_rows(COG20_occurrence,KOfam_occurrence)
  
  # Initialize an empty data frame for the final table
  final_table <- data.frame()
  
  # Iterate over each gene in Gene_Info
  for (i in 1:nrow(Gene_Info)) {
    gene_name <- Gene_Info$Gene_name[i]
    kofam_acc <- Gene_Info$KOfam_acc[i]
    cog_acc <- Gene_Info$COG_acc[i]
    
    # Filter occurrences for the current gene based on KOfam_acc or COG_acc
    gene_occurrences <- occurrences %>%
      filter(accession == kofam_acc | accession == cog_acc) %>%
      select(genome_name, gene_callers_id, gene_function)
    
    # Add the gene name to the occurrences
    gene_occurrences <- gene_occurrences %>%
      mutate(gene = gene_name)
    
    # Merge with Genomes to get the Group_ID
    gene_occurrences <- gene_occurrences %>%
      inner_join(Genomes, by = c("genome_name" = "Pangenome_ID")) %>%
      select(genome_name, Group_ID, gene, gene_callers_id)
    
    # Bind to the final table
    final_table <- bind_rows(final_table, gene_occurrences)
  }
  
  # Remove duplicate gene and gene_callers_id entries
  final_table <- final_table %>%
    distinct(genome_name, Group_ID, gene, gene_callers_id, .keep_all = TRUE)
  
  return(final_table)
}

# Example usage
final_table <- generate_gene_callers_list(COG20_occurrence, KOfam_occurrence, Genomes, Gene_Info)


# Merge results with genome habitat preference results
strain_level_habitat_classification_res <- read.table("/Users/home/SPECIES_LEVEL_PANGENOMES/TESTS/Classifying_genomes_to_oral_sites/Groups/Without_HP/Final_results_strain_level.txt", header = TRUE, sep = "\t")

merged_res <- merge(final_table, strain_level_habitat_classification_res, by.x = "genome_name", by.y = "Pangenome_ID")

# clean up
merged_res <- merged_res %>% 
  dplyr::select(-c(Group_ID.y)) %>% 
  dplyr::rename(Group_ID = Group_ID.x) 

# Save the occurrence table
write.table(merged_res, "/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/Nitrate/gene_callers_list.txt", sep = "\t", row.names = FALSE, quote = FALSE)


```


Send gene caller ID metadata to Barhal...
```{bash, eval = FALSE}

scp -r /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/Nitrate/gene_callers_list.txt BARHAL_SERVER:/workspace/jmarkwelchlab/P_0003_Neisseriaceae/DENITRIFICATION_genes/gene_callers_list.txt

```

Get ECG status for each denitrification gene for the three main sites.
```{bash, eval = FALSE}

# Directories and file paths
MAIN_DIR="/workspace/jmarkwelchlab/P_0003_Neisseriaceae"
GENE_DIR="$MAIN_DIR/22_GENE_LEVEL"
metadata="$MAIN_DIR/DENITRIFICATION_genes/gene_callers_list.txt"
OUT="$MAIN_DIR/DENITRIFICATION_genes/DENITRIFICATION_genes_ECG_status.txt"


# Create the output file and write the header
echo -e "gene_caller_id\tTD\tSUPP\tKG" > $OUT

# Read the metadata file and process each geneome
while read genome_name gene_caller_id; do
  # Initialize detection status for each site as "NA" (or any placeholder for missing data)
  TD="NA"
  SUPP="NA"
  KG="NA"

  # Loop through each site and get the detection status
  for site in TD SUPP KG; do
    # Check if the detection file exists
    DETECTION_FILE="$GENE_DIR/$site/${genome_name}-ENV-DETECTION.txt"
    if [[ -f $DETECTION_FILE ]]; then
      # Get detection status for the current gene caller ID
      DETECTION=$(grep -w "$gene_caller_id" "$DETECTION_FILE" | awk '{print $2}')
      # Update the detection status variable for the current site
      if [[ $site == "TD" ]]; then
        TD=$DETECTION
      elif [[ $site == "SUPP" ]]; then
        SUPP=$DETECTION
      elif [[ $site == "KG" ]]; then
        KG=$DETECTION
      fi
    fi
  done

  # Write the results to the output file
  echo -e "$gene_caller_id\t$TD\t$SUPP\t$KG" >> $OUT

done < <(awk 'NR>1 {print $1, $4}' $metadata)



```

Merge ECG status file with Hia gene metadata

```{r, eval = FALSE}

metadata <- read.table("/workspace/jmarkwelchlab/P_0003_Neisseriaceae/DENITRIFICATION_genes/gene_callers_list.txt", header = TRUE, sep = "\t")

ECG_status <- read.table("/workspace/jmarkwelchlab/P_0003_Neisseriaceae/DENITRIFICATION_genes/DENITRIFICATION_genes_ECG_status.txt", header = TRUE, sep = "\t")

merged_file <- merge(metadata, ECG_status, by.x = "gene_callers_id", by.y = "gene_caller_id")

write.table(merged_file,"/workspace/jmarkwelchlab/P_0003_Neisseriaceae/DENITRIFICATION_genes/DENITRIFICATION_genes_metadata_with_ECG.txt", row.names = FALSE, quote = FALSE, sep = "\t")
```

Send results to local...
```{bash, eval = FALSE}

scp -r BARHAL_SERVER:/workspace/jmarkwelchlab/P_0003_Neisseriaceae/DENITRIFICATION_genes/DENITRIFICATION_genes_metadata_with_ECG.txt /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/Nitrate/DENITRIFICATION_genes_metadata_with_ECG.txt

```

Get mean and median coverage across samples for each denitrification gene for the three main sites.

```{r, eval = FALSE}

library("dplyr")
library(readr)

# Function to calculate mean, standard deviation, and median
calculate_stats <- function(values) {
  mean_val <- mean(values, na.rm = TRUE)
  sd_mean_val <- sd(values, na.rm = TRUE)
  median_val <- median(values, na.rm = TRUE)
  
  return(c(mean_val, sd_mean_val, median_val))
}

# Main directory
MAIN_DIR <- "/workspace/jmarkwelchlab/P_0003_Neisseriaceae"
GENE_DIR <- file.path(MAIN_DIR, "22_GENE_LEVEL")
metadata <- file.path(MAIN_DIR, "DENITRIFICATION_genes/gene_callers_list.txt")
OUT <- file.path(MAIN_DIR, "DENITRIFICATION_genes/DENITRIFICATION_genes_coverage_data.txt")

# Create the output file and write the header
write("gene_callers_id\tTD_mean\tTD_sd_mean\tTD_median\tSUPP_mean\tSUPP_sd_mean\tSUPP_median\tKG_mean\tKG_sd_mean\tKG_median", file = OUT)

# Read the metadata file
metadata_df <- read_delim(metadata, delim = "\t", col_names = TRUE)


# Process each genome
for (i in 1:nrow(metadata_df)) {
  genome_name <- metadata_df[i, 1]
  gene_callers_id <- metadata_df[i, 4]
  
  # Initialize detection status for each site as "NA"
  TD_mean <- TD_sd_mean <- TD_median <- "NA"
  SUPP_mean <- SUPP_sd_mean <- SUPP_median <- "NA"
  KG_mean <- KG_sd_mean <- KG_median <- "NA"
  
  # Loop through each site and get the coverage
  for (site in c("TD", "SUPP", "KG")) {
    COVERAGE_FILE <- file.path(GENE_DIR, site, paste0(genome_name, "-GENE-COVs.txt"))
    
    if (file.exists(COVERAGE_FILE)) {
      # Read the coverage file
      coverage_data <- read_delim(COVERAGE_FILE, delim = "\t", col_names = TRUE)
      
      # Get all coverage values for the current gene caller ID
      coverage_values <- coverage_data %>% filter(key == gene_callers_id$gene_callers_id) %>% select(-key) %>% unlist()
      
      # Check if coverage values are not empty
      if (length(coverage_values) > 0) {
        # Calculate statistics
        stats <- calculate_stats(as.numeric(coverage_values))
        
        # Update the corresponding variables
        if (site == "TD") {
          TD_mean <- stats[1]
          TD_sd_mean <- stats[2]
          TD_median <- stats[3]
        } else if (site == "SUPP") {
          SUPP_mean <- stats[1]
          SUPP_sd_mean <- stats[2]
          SUPP_median <- stats[3]
        } else if (site == "KG") {
          KG_mean <- stats[1]
          KG_sd_mean <- stats[2]
          KG_median <- stats[3]
        }
      }
    }
  }
  
  # Write the results to the output file
  result <- paste(gene_callers_id, TD_mean, TD_sd_mean, TD_median,
                  SUPP_mean, SUPP_sd_mean, SUPP_median,
                  KG_mean, KG_sd_mean, KG_median, sep = "\t")
  write(result, file = OUT, append = TRUE)
}

```

Merge gene coverage file with gene metadata

```{r, eval = FALSE}

metadata <- read.table("/workspace/jmarkwelchlab/P_0003_Neisseriaceae/DENITRIFICATION_genes/gene_callers_list.txt", header = TRUE, sep = "\t")

coverage <- read.table("/workspace/jmarkwelchlab/P_0003_Neisseriaceae/DENITRIFICATION_genes/DENITRIFICATION_genes_coverage_data.txt", header = TRUE, sep = "\t")

merged_file <- merge(metadata, coverage, by = "gene_callers_id")

write.table(merged_file,"/workspace/jmarkwelchlab/P_0003_Neisseriaceae/DENITRIFICATION_genes/DENITRIFICATION_genes_metadata_with_coverage.txt", row.names = FALSE, quote = FALSE, sep = "\t")

# Check for duplicates
duplicate_gene_callers_id <- merged_file %>%
  group_by(gene_callers_id) %>%
  filter(n() > 1) %>%
  ungroup()

# View the rows with duplicate gene_callers_id values
print(duplicate_gene_callers_id)
```


Send results to local...
```{bash, eval = FALSE}

scp -r BARHAL_SERVER:/workspace/jmarkwelchlab/P_0003_Neisseriaceae/DENITRIFICATION_genes/DENITRIFICATION_genes_metadata_with_coverage.txt /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/Nitrate/DENITRIFICATION_genes_metadata_with_coverage.txt

```

Now I want to merge the detection and coverage data with the gene occurance data... We will start by doing so for a subset of the genes, those required to encode nitrate reductase, narGHIJ.

```{r, eval = FALSE}

occurrence_table <- read_delim("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/Nitrate/DENITRIFICATION_genes_metadata_with_ECG.txt", delim = "\t")

ecg_table <- read_delim("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/Nitrate/DENITRIFICATION_genes_metadata_with_ECG.txt" , delim = "\t")

coverage_table <- read_delim("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/Nitrate/DENITRIFICATION_genes_metadata_with_coverage.txt", delim = "\t")


# Remove duplicate rows in the ecg_table and coverage_table, keeping the first occurrence
ecg_table <- ecg_table %>% distinct(gene_callers_id, .keep_all = TRUE)
coverage_table <- coverage_table %>% distinct(gene_callers_id, .keep_all = TRUE)
  


merge_files <- function(occurrence_file, ecg_file, coverage_file, output_file) {
  # Load necessary libraries
  library("dplyr")
  library(readr)
  
  # Read the input files
  occurrence_table <- read_delim(occurrence_file, delim = "\t")
  ecg_table <- read_delim(ecg_file, delim = "\t")
  coverage_table <- read_delim(coverage_file, delim = "\t")
  
  # Remove duplicate rows in the ecg_table and coverage_table, keeping the first occurrence
  ecg_table <- ecg_table %>% distinct(gene_callers_id, .keep_all = TRUE)
  coverage_table <- coverage_table %>% distinct(gene_callers_id, .keep_all = TRUE)
  
  # Select necessary columns from the ECG table for merging
  ecg_selected <- ecg_table %>%
    select(gene_callers_id, gene, genome_name, Group_ID, TD, SUPP, KG) %>%
    rename(TD_DETECTION = TD, SUPP_DETECTION = SUPP, KG_DETECTION = KG)
  
  # Select necessary columns from the coverage table for merging
  coverage_selected <- coverage_table %>%
    select(gene_callers_id, genome_name, TD_mean, TD_sd_mean, TD_median,
           SUPP_mean, SUPP_sd_mean, SUPP_median,
           KG_mean, KG_sd_mean, KG_median)
  
  # Initialize the merged data frame with occurrence_table columns
  merged_df <- occurrence_table %>%
    select(Pangenome_ID, Group_ID, ABUNDANT_IN, PREVALENT_IN, COMBINED_CLASSIFICATION, GENERAL_SITE_CLASS)
  
  # List of target genes
  target_genes <- c("narG", "narH", "narI", "narJ")
  
  for (gene in target_genes) {
    # Merge with ECG table
    ecg_gene <- ecg_selected %>%
      filter(gene == !!gene) %>%
      select(genome_name, gene_callers_id, Group_ID, TD_DETECTION, SUPP_DETECTION, KG_DETECTION) %>% 
      rename(Pangenome_ID = genome_name)
    
    coverage_gene <- coverage_selected %>%
      filter(gene_callers_id %in% ecg_selected$gene_callers_id[ecg_selected$gene == gene]) %>%
      rename(Pangenome_ID = genome_name) %>%
      select(-gene_callers_id)
    
    # Merge with occurrence table
    occurrence_gene <- occurrence_table %>%
      select(Pangenome_ID, Group_ID, gene) %>%
      left_join(ecg_gene, by = c("Pangenome_ID","Group_ID")) %>%
      left_join(coverage_gene, by = "Pangenome_ID")
    
    # Rename columns for merging
    colnames(occurrence_gene)[3] <- gene
    colnames(occurrence_gene)[4] <- paste0(gene,"_gene_caller_id")
    colnames(occurrence_gene)[5:7] <- paste0(gene, "_", c("TD_DETECTION", "SUPP_DETECTION", "KG_DETECTION"))
    colnames(occurrence_gene)[8:16] <- paste0(gene, "_", c("TD_mean", "TD_sd_mean", "TD_median",
                                                     "SUPP_mean", "SUPP_sd_mean", "SUPP_median",
                                                     "KG_mean", "KG_sd_mean", "KG_median"))
    
    # Merge into the main dataframe
    merged_df <- merged_df %>%
      left_join(occurrence_gene, by = c("Pangenome_ID", "Group_ID"))
  }
  
  # Write the output file
  write_delim(merged_df, output_file, delim = "\t")
}

# Usage
merge_files("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/Nitrate/occurrence_table.txt", "/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/Nitrate/DENITRIFICATION_genes_metadata_with_ECG.txt", "/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/Nitrate/DENITRIFICATION_genes_metadata_with_coverage.txt",
"/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/Nitrate/Nitrate_reductase_occurrence_and_ECG_coverage.txt")

merge_files_without_coverage <- function(occurrence_file, ecg_file, output_file) {
  # Load necessary libraries
  library("dplyr")
  library(readr)
  
  # Read the input files
  occurrence_table <- read_delim(occurrence_file, delim = "\t")
  ecg_table <- read_delim(ecg_file, delim = "\t")
  
  # Remove duplicate rows in the ecg_table, keeping the first occurrence
  ecg_table <- ecg_table %>% distinct(gene_callers_id, genome_name, gene, .keep_all = TRUE)
  
  # Select necessary columns from the ECG table for merging
  ecg_selected <- ecg_table %>%
    select(gene_callers_id, gene, genome_name, Group_ID, TD, SUPP, KG) %>%
    rename(TD_DETECTION = TD, SUPP_DETECTION = SUPP, KG_DETECTION = KG)
  
  # Initialize the merged data frame with occurrence_table columns
  merged_df <- occurrence_table %>%
    select(Pangenome_ID, Group_ID, ABUNDANT_IN, PREVALENT_IN, COMBINED_CLASSIFICATION, GENERAL_SITE_CLASS)
  
  # List of target genes
  target_genes <- c("narG", "narH", "narI", "narJ")
  
  for (gene in target_genes) {
    # Merge with ECG table
    ecg_gene <- ecg_selected %>%
      filter(gene == !!gene) %>%
      select(gene_callers_id, genome_name, Group_ID, TD_DETECTION, SUPP_DETECTION, KG_DETECTION) %>% 
      rename(Pangenome_ID = genome_name)
    
    # Merge with occurrence table
    occurrence_gene <- occurrence_table %>%
      select(Pangenome_ID, Group_ID, gene)
    
    gene_data <- occurrence_gene %>%
      left_join(ecg_gene, by = c("Pangenome_ID", "Group_ID"))
    
    # Rename columns for merging
    colnames(gene_data)[3] <- gene
    colnames(gene_data)[4] <- paste0(gene,"_gene_caller_id")
    colnames(gene_data)[5:7] <- paste0(gene, "_", c("TD_DETECTION", "SUPP_DETECTION", "KG_DETECTION"))
    
    # Merge into the main data frame
    merged_df <- merged_df %>%
      left_join(gene_data, by = c("Pangenome_ID", "Group_ID"))
  }
  
  # Write the output file
  write_delim(merged_df, output_file, delim = "\t")
}

# Usage
merge_files_without_coverage("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/Nitrate/occurrence_table.txt", "/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/Nitrate/DENITRIFICATION_genes_metadata_with_ECG.txt",  "/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/Nitrate/Nitrate_reductase_occurrence_and_ECG.txt")

```




# 17 Functional occurence plots

## 17.1 COG20 Functional occurrence plot


```{r, eval = FALSE}


library("tibble")
library("dplyr")
library("tidyr")
library("pheatmap")
library("dendextend")


# load data
COG20_occurence_df <- read.table("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/20_FUNCTIONAL_ENRICHMENT/COG20_FUNCTION_Pangenome_FUNC_OCCURRENCE_GRANULAR.TXT", sep = "\t", header = TRUE, fill = TRUE, quote = "")

# change column name for first column
COG20_occurence_df <- COG20_occurence_df %>% 
  rename(COG20_function = X)

# load genome metadata 
habitat_prefs <- read.table("/Users/home/SPECIES_LEVEL_PANGENOMES/TESTS/Classifying_genomes_to_oral_sites/Groups/Without_HP/Refined_habitat_prefs.txt", sep = "\t", header = TRUE, fill = TRUE, quote = "", stringsAsFactors = TRUE)  # Make sure to set the correct separator if it's not tab.

habitat_prefs <- habitat_prefs %>% 
  mutate(Genus = case_when(grepl('N_', Group_ID) ~ 'Neisseria',
                           grepl('K_', Group_ID) ~ 'Kingella',
                           grepl('E_', Group_ID) ~ 'Eikenella',
                           grepl('S_', Group_ID) ~ 'Simonsiella')) %>% 
  mutate(Species = case_when(grepl('E_corrodens', Group_ID) ~ 'corrodens',
                             grepl('E_exigua', Group_ID) ~ 'exigua',
                             grepl('E_glucosivorans', Group_ID) ~ 'glucosivorans',
                             grepl('E_halliae', Group_ID) ~ 'halliae',
                             grepl('K_bonacorsii_oralis', Group_ID) ~ 'oralis',
                             grepl('K_denitrificans', Group_ID) ~ 'denitrificans',
                             grepl('K_kingae', Group_ID) ~ 'kingae',
                             grepl('K_negevensis', Group_ID) ~ 'negevensis',
                             grepl('K_potus', Group_ID) ~ 'potus',
                             grepl('N_maigaei', Group_ID) ~ 'maigaei',
                             grepl('K_sp_str_SNUBH_2017', Group_ID) ~ 'sp_str_SNUBH_2017',
                             grepl('N_bacilliformis', Group_ID) ~ 'bacilliformis',
                             grepl('N_basseii', Group_ID) ~ 'basseii',
                             grepl('N_benedictiae', Group_ID) ~ 'benedictiae',
                             grepl('N_bergeri', Group_ID) ~ 'bergeri',
                             grepl('N_blantyrii', Group_ID) ~ 'blantyrii',
                             grepl('N_brasiliensis', Group_ID) ~ 'brasiliensis',
                             grepl('N_cinerea', Group_ID) ~ 'cinerea',
                             grepl('N_cinerea_str_CCUG_5746', Group_ID) ~ 'cinerea_str_CCUG_5746',
                             grepl('N_dentiae', Group_ID) ~ 'dentiae',
                             grepl('N_dumasiana', Group_ID) ~ 'dumasiana',
                             grepl('N_elongata', Group_ID) ~ 'elongata',
                             grepl('N_gonorrhoeae', Group_ID) ~ 'gonorrhoeae',
                             grepl('N_lactamica', Group_ID) ~ 'lactamica',
                             grepl('N_lactamica_str_NS19', Group_ID) ~ 'lactamica_str_NS19',
                             grepl('N_mucosa_', Group_ID) ~ 'mucosa',
                             grepl('N_oralis', Group_ID) ~ 'oralis',
                             grepl('N_polysaccharea', Group_ID) ~ 'polysaccharea',
                             grepl('N_shayeganii', Group_ID) ~ 'shayeganii',
                             grepl('N_sp_HMT_020', Group_ID) ~ 'sp_HMT_020',
                             grepl('N_sp_str_3986_and_str_51_81', Group_ID) ~ 'sp_str_3986_and_str_51_81',
                             grepl('N_sp_str_HSC_16F19', Group_ID) ~ 'sp_str_HSC_16F19',
                             grepl('N_sp_str_MVDL20_010259', Group_ID) ~ 'sp_str_MVDL20_010259',
                             grepl('N_subflava_I', Group_ID) ~ 'subflava_A',
                             grepl('N_subflava_II', Group_ID) ~ 'subflava_A',
                             grepl('N_subflava_III', Group_ID) ~ 'subflava_A',
                             grepl('N_subflava_IV', Group_ID) ~ 'subflava_A',
                             grepl('N_subflava_VII', Group_ID) ~ 'subflava_A',
                             grepl('N_subflava_VIII', Group_ID) ~ 'subflava_A',
                             grepl('N_subflava_V', Group_ID) ~ 'subflava_B',
                             grepl('N_uirgultaei', Group_ID) ~ 'uirgultaei',
                             grepl('N_viridiae', Group_ID) ~ 'viridiae',
                             grepl('N_wadsworthii', Group_ID) ~ 'wadsworthii',
                             grepl('N_weaveri', Group_ID) ~ 'weaveri',
                             grepl('N_meningitidis', Group_ID) ~ 'meningitidis',
                             grepl('S_muelleri', Group_ID) ~ 'muelleri'))




##### Remove undetected genomes from main df, which is any pangenome ID in the habitat_prefs df that has a granular preference of "NONE"
#make list of detected genomes
detected_genomes <- habitat_prefs %>% 
  filter(!(Habitat_Preference_Granular == "NONE")) %>% 
  select(item)
# filter columns
COG20_occurence_df_filtered <- COG20_occurence_df %>% 
  select(COG20_function, any_of(detected_genomes$item))


##### 
COG20_occurence_df_filtered_binary <- COG20_occurence_df_filtered %>%
  mutate(across(-COG20_function, ~ ifelse(. >= 1, 1, 0))) %>% 
  mutate(COG20_function_seq = seq(1, length(COG20_occurence_df_filtered$COG20_function))) %>% 
  select(COG20_function, COG20_function_seq, everything())

# save functions and seqs
Functions_and_numbers <- COG20_occurence_df_filtered_binary %>% 
  select(COG20_function, COG20_function_seq)


#### Make final matrix for clustering
COG20_occurence_df_filtered_binary_final <- COG20_occurence_df_filtered_binary %>% 
  select(everything(), -COG20_function )

# rename row names to species IDs
COG20_occurence_df_filtered_binary_final_2 <- COG20_occurence_df_filtered_binary_final[,-1]
rownames(COG20_occurence_df_filtered_binary_final_2) <- COG20_occurence_df_filtered_binary_final[,1]

# convert to matrix
COG20_occurence_matrix <- as.matrix(COG20_occurence_df_filtered_binary_final_2)

# Identify empty rows (rows that only contain zeros)
empty_rows <- apply(COG20_occurence_matrix, 1, function(x) all(x == 0))

# Optionally, print which rows are empty
print(which(empty_rows))

# Remove empty rows from the matrix
COG20_occurence_matrix_no_empty <- COG20_occurence_matrix[!empty_rows, , drop = FALSE]

# transmutate so genomes are Row ids and cog functions are columns
trans_COG20_occurence_matrix_no_empty <- t(COG20_occurence_matrix_no_empty)

# Create a mapping for multiple annotations
annotation_row <- habitat_prefs %>%
  select(item, Habitat_Preference_Granular, Species, Genus) %>%
  rename(Genome = item, 'Habitat Preference' = Habitat_Preference_Granular) %>%
  filter(Genome %in% rownames(trans_COG20_occurence_matrix_no_empty)) %>%
  column_to_rownames(var = "Genome")  # Set genome IDs as rownames

# Specify colors
ann_colors = list(
  'Habitat Preference' = c(PLAQUE = "blue", TONGUE= "red", GUMS = "green"),
  Species = c(bacilliformis = "#c41c00", cinerea = "#06ea16", corrodens = "#01278f",denitrificans = "#0020c2",
              elongata = "#70eafa",exigua = "#f7e702",halliae = "#c2045d",mucosa = "#914747",
              muelleri = "#9a5bb5",oralis = "#ff00dd",sp_HMT_020 = "#8200f4",sp_str_SNUBH_2017 = "#f77cd2",
              subflava = "#ff006f"),
  Genus = c(Eikenella = "#6363ea", Kingella = "#25f4b9", Neisseria = "#f46c6c", Simonsiella = "#edcb3b")
)


save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# Assuming 'trans_COG20_occurence_matrix_no_empty' is your data matrix
# and 'annotation_row' is your annotation data frame

# Step 1: Create initial pheatmap to extract clustering information
initial_plot <- pheatmap(
  trans_COG20_occurence_matrix_no_empty,
  cluster_rows = TRUE,      
  cluster_cols = TRUE, 
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  color = c("black", "cyan"), # Specify colors for binary heatmap: 0 = white, 1 = black
  show_rownames = FALSE,     # Hide row names if desired
  show_colnames = FALSE,
  annotation_row = annotation_row, # Hide column names if desired
  annotation_colors = ann_colors,
  legend = FALSE,
  border_color = "black",
  annotation_legend = FALSE,
  treeheight_col = 100
)

# Step 2: Extract row and column dendrograms
row_tree <- as.dendrogram(initial_plot$tree_row)
col_tree <- as.dendrogram(initial_plot$tree_col)

# Step 3: Ladderize the dendrograms
ladderized_row_tree <- ladderize(row_tree)
ladderized_col_tree <- ladderize(col_tree)

# Step 4: Extract the new order from ladderized trees
row_order <- order.dendrogram(ladderized_row_tree)
#col_order <- order.dendrogram(ladderized_col_tree)

ladderized_col_phylo_tree <- as.phylo(ladderized_col_tree)
ggtree_ladderized_col_tree <- ggtree(ladderized_col_phylo_tree, ladderize=F) +
  geom_nodepoint() +
  geom_text(aes(label=node), hjust=-.3)
rotated_ggtree_ladderized_col_tree <- rotate(ggtree_ladderized_col_tree, 2297)
tip_labels <- rev(get_taxa_name(rotated_ggtree_ladderized_col_tree))


# Step 5: Reorder the data matrix
data_matrix_ladderized <- trans_COG20_occurence_matrix_no_empty[row_order, tip_labels]

# Step 6: Plot pheatmap with ladderized order
ladderized_plot <- pheatmap(
  data_matrix_ladderized,
  cluster_rows = FALSE,       # Rows are already ladderized
  cluster_cols = FALSE,       # Columns are already ladderized
  color = c("black", "cyan"), # Specify colors for binary heatmap: 0 = white, 1 = cyan
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_row = annotation_row, # Hide column names if desired
  annotation_colors = ann_colors,
  annotation_names_row = FALSE,
  legend = FALSE,
  border_color = NA,
  annotation_legend = FALSE
)

save_pheatmap_pdf(ladderized_plot, "/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/FIGURES/Functional_enrichment/ladderized_plot.pdf", width = 11, height = 8)


############################################################################################
############################################################################################

library("tidyr")
library("dplyr")

# Extract relevant information from habitat_prefs
genome_metadata <- habitat_prefs %>%
  filter(item %in% colnames(COG20_occurence_df_filtered)) %>% 
  select(item, Genus, Species, Habitat_Preference_Granular) 

# Convert COG20_occurence_df_filtered to long format
COG20_occurence_long <- COG20_occurence_df_filtered %>%
  pivot_longer(
    cols = -COG20_function,
    names_to = "Genome",
    values_to = "Presence"
  )

# Join metadata to add genus-species information
COG20_occurence_long <- COG20_occurence_long %>%
  left_join(genome_metadata, by = c("Genome" = "item"))

# Collapse data to genus-species level
COG20_occurence_habitat <- COG20_occurence_long %>%
  group_by(COG20_function, Habitat_Preference_Granular) %>%
  summarize(Presence = sum(Presence, na.rm = TRUE)) %>%  # Use max to indicate presence/absence at genus-species level
  pivot_wider(
    names_from = Habitat_Preference_Granular,
    values_from = Presence,
    values_fill = list(Presence = 0)
  )


COG20_occurence_habitat <- as.data.frame(COG20_occurence_habitat)

COG20_occurence_habitat_binary <- COG20_occurence_habitat %>%
  mutate(across(-COG20_function, ~ ifelse(. >= 1, 1, 0))) %>% 
  mutate(COG20_function_seq = seq(1, length(COG20_occurence_habitat$COG20_function))) %>% 
  select(COG20_function, COG20_function_seq, everything())

#### Make final matrix for clustering
COG20_occurence_habitat_binary_final <- COG20_occurence_habitat_binary %>% 
  select(everything(), -COG20_function )

# rename row names to species IDs
COG20_occurence_habitat_binary_final_2 <- COG20_occurence_habitat_binary_final[,-1]
rownames(COG20_occurence_habitat_binary_final_2) <- COG20_occurence_habitat_binary_final[,1]

# convert to matrix
COG20_occurence_habitat_matrix <- as.matrix(COG20_occurence_habitat_binary_final_2)

# Identify empty rows (rows that only contain zeros)
COG20_occurence_habitat_empty_rows <- apply(COG20_occurence_habitat_matrix, 1, function(x) all(x == 0))

# Remove empty rows from the matrix
COG20_occurence_habitat_matrix_no_empty <- COG20_occurence_habitat_matrix[!COG20_occurence_habitat_empty_rows, , drop = FALSE]

# transmutate so genomes are Row ids and cog functions are columns
trans_COG20_occurence_habitat_matrix_no_empty <- t(COG20_occurence_habitat_matrix_no_empty)




# Specify colors
habitat_ann_colors = list(
  'Habitat Preference' = c(PLAQUE = "blue", TONGUE= "red", GUMS = "green"),
  Species = c(X = "white"),
  Genus = c(X = "white")
)

habitat_annotation_row <- genome_metadata %>%
  select(Habitat_Preference_Granular) %>%
  mutate('Habitat Preference' = Habitat_Preference_Granular,
         Species = "X", 
         Genus = "X") %>%
  distinct() %>% 
  column_to_rownames(var = "Habitat_Preference_Granular")  # Set genome IDs as rownames


# Step 1: Create initial pheatmap to extract clustering information
habitat_initial_plot <- pheatmap(
  trans_COG20_occurence_habitat_matrix_no_empty,
  cluster_rows = TRUE,      
  cluster_cols = TRUE, 
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  color = c("black", "cyan"),
  show_rownames = TRUE,     # Hide row names if desired
  show_colnames = FALSE,
  annotation_row = habitat_annotation_row, # Hide column names if desired
  annotation_colors = habitat_ann_colors,
  annotation_names_row = FALSE,
  legend = FALSE,
  border_color = "black",
  annotation_legend = FALSE
)

# Step 2: Extract row and column dendrograms
habitat_row_tree <- as.dendrogram(habitat_initial_plot$tree_row)
habitat_col_tree <- as.dendrogram(habitat_initial_plot$tree_col)

# Step 3: Ladderize the dendrograms
habitat_ladderized_row_tree <- ladderize(habitat_row_tree)
habitat_ladderized_col_tree <- ladderize(habitat_col_tree)

# Step 4: Extract the new order from ladderized trees
habitat_row_order <- rev(order.dendrogram(habitat_ladderized_row_tree))
#col_order <- order.dendrogram(ladderized_col_tree)

habitat_ladderized_col_phylo_tree <- as.phylo(habitat_ladderized_col_tree)
habitat_ggtree_ladderized_col_tree <- ggtree(habitat_ladderized_col_phylo_tree, ladderize=F) +
  geom_nodepoint() +
  geom_text(aes(label=node), hjust=-.3)
habitat_rotated_ggtree_ladderized_col_tree <- rotate(habitat_ggtree_ladderized_col_tree, 2297)
habitat_tip_labels <- rev(get_taxa_name(habitat_rotated_ggtree_ladderized_col_tree))


# Step 5: Reorder the data matrix
habitat_data_matrix_ladderized <- trans_COG20_occurence_habitat_matrix_no_empty[habitat_row_order, habitat_tip_labels]

# Step 6: Plot pheatmap with ladderized order
habitat_ladderized_plot <- pheatmap(
  habitat_data_matrix_ladderized,
  cluster_rows = FALSE,       # Rows are already ladderized
  cluster_cols = FALSE,       # Columns are already ladderized
  color = c("black", "cyan"), # Specify colors for binary heatmap: 0 = white, 1 = cyan
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_row = habitat_annotation_row, # Hide column names if desired
  annotation_colors = habitat_ann_colors,
  annotation_names_row = FALSE,
  legend = FALSE,
  border_color = NA,
  annotation_legend = FALSE
)

save_pheatmap_pdf(habitat_ladderized_plot, "/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/FIGURES/Functional_enrichment/habitat_ladderized_plot.pdf", width = 11, height = 4)


############################################################################################
############################################################################################

library(ggplotify)
habitat_ladderized_plot_ggplot <- as.ggplot(habitat_ladderized_plot)
ladderized_plot_ggplot <- as.ggplot(ladderized_plot)

combo <- egg::ggarrange(ladderized_plot_ggplot, habitat_ladderized_plot_ggplot, ncol = 1, heights = c(1,0.2))
ggsave("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/FIGURES/Functional_enrichment/COMBO_ladderized_plot.pdf", plot = combo, height = 8, width = 10)


############################################################################################
############################################################################################

ladderized_row_phylo_tree <- as.phylo(ladderized_row_tree)
ggtree_ladderized_row_tree <- ggtree(ladderized_row_phylo_tree, ladderize=F)

ggsave("/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/FIGURES/Functional_enrichment/ggtree_ladderized_row_tree.pdf", plot = ggtree_ladderized_row_tree, height = 6.4, width = 2)





```



Species-level COG20 fucntional annotations binary map...

```{r, eval = FALSE}


library("tidyr")
library("dplyr")

# Extract relevant information from habitat_prefs
genome_metadata <- habitat_prefs %>%
  filter(item %in% colnames(COG20_occurence_df_filtered)) %>% 
  select(item, Genus, Species, Habitat_Preference_Granular) %>% 
  mutate(Genus_Species = paste(Genus, Species, sep = "_"))  # Create a combined genus-species identifier


# Convert COG20_occurence_df_filtered to long format
COG20_occurence_long <- COG20_occurence_df_filtered %>%
  pivot_longer(
    cols = -COG20_function,
    names_to = "Genome",
    values_to = "Presence"
  )

# Join metadata to add genus-species information
COG20_occurence_long <- COG20_occurence_long %>%
  left_join(genome_metadata, by = c("Genome" = "item"))

# Collapse data to genus-species level
COG20_occurence_collapsed <- COG20_occurence_long %>%
  group_by(COG20_function, Genus_Species) %>%
  summarize(Presence = sum(Presence, na.rm = TRUE)) %>%  # Use max to indicate presence/absence at genus-species level
  pivot_wider(
    names_from = Genus_Species,
    values_from = Presence,
    values_fill = list(Presence = 0)
  )


COG20_occurence_collapsed <- as.data.frame(COG20_occurence_collapsed)

COG20_occurence_collapsed_binary <- COG20_occurence_collapsed %>%
  mutate(across(-COG20_function, ~ ifelse(. >= 1, 1, 0))) %>% 
  mutate(COG20_function_seq = seq(1, length(COG20_occurence_collapsed$COG20_function))) %>% 
  select(COG20_function, COG20_function_seq, everything())

#### Make final matrix for clustering
COG20_occurence_collapsed_binary_final <- COG20_occurence_collapsed_binary %>% 
  select(everything(), -COG20_function )

# rename row names to species IDs
COG20_occurence_collapsed_binary_final_2 <- COG20_occurence_collapsed_binary_final[,-1]
rownames(COG20_occurence_collapsed_binary_final_2) <- COG20_occurence_collapsed_binary_final[,1]

# convert to matrix
COG20_occurence_collapsed_matrix <- as.matrix(COG20_occurence_collapsed_binary_final_2)

# Identify empty rows (rows that only contain zeros)
COG20_occurence_collapsed_empty_rows <- apply(COG20_occurence_collapsed_matrix, 1, function(x) all(x == 0))

# Remove empty rows from the matrix
COG20_occurence_collapsed_matrix_no_empty <- COG20_occurence_collapsed_matrix[!COG20_occurence_collapsed_empty_rows, , drop = FALSE]

# transmutate so genomes are Row ids and cog functions are columns
trans_COG20_occurence_collapsed_matrix_no_empty <- t(COG20_occurence_collapsed_matrix_no_empty)




# Create a mapping for multiple annotations
species_annotation_row <- genome_metadata %>%
  select(Genus_Species, Habitat_Preference_Granular, Species, Genus) %>%
  rename('Habitat Preference' = Habitat_Preference_Granular) %>%
  filter(Genus_Species %in% rownames(trans_COG20_occurence_collapsed_matrix_no_empty)) %>%
  distinct() %>% 
  column_to_rownames(var = "Genus_Species")  # Set genome IDs as rownames

# Specify colors
species_ann_colors = list(
  'Habitat Preference' = c(PLAQUE = "blue", TONGUE= "red", GUMS = "green"),
  Species = c(bacilliformis = "#c41c00", cinerea = "#06ea16", corrodens = "#01278f",denitrificans = "#0020c2",
              elongata = "#70eafa",exigua = "#f7e702",halliae = "#c2045d",mucosa = "#914747",
              muelleri = "#9a5bb5",oralis = "#ff00dd",sp_HMT_020 = "#8200f4",sp_str_SNUBH_2017 = "#f77cd2",
              subflava_A = "#ff006f", subflava_B = "orange"),
  Genus = c(Eikenella = "#6363ea", Kingella = "#25f4b9", Neisseria = "#f46c6c", Simonsiella = "#edcb3b")
)

species_occurrence_plot <- pheatmap(
  trans_COG20_occurence_collapsed_matrix_no_empty,
  cluster_rows = TRUE,      
  cluster_cols = TRUE, 
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  color = c("black", "cyan"),
  show_rownames = TRUE,     # Hide row names if desired
  show_colnames = FALSE,
  annotation_row = species_annotation_row, # Hide column names if desired
  annotation_colors = species_ann_colors,
  legend = FALSE,
  border_color = "black",
  annotation_legend = FALSE,
  treeheight_col = 100
)

save_pheatmap_pdf(species_occurrence_plot, "/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/FIGURES/Functional_enrichment/test_3.pdf", width = 11, height = 8)




```



## 17.2 KEGG Modules occurrence plot

```{r, eval = FALSE}


save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

library("tibble")
library("dplyr")
library("tidyr")
library("pheatmap")
library("dendextend")



MainDIR <- "/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/"

# load data
MetabolicPaths_modules <- read.table(paste0(MainDIR,"12_FUNCTIONAL_ANNOTATION/KEGG_Metabolic_pathways/P_0003_Neisseriaceae-MetabolicPaths_modules.txt"), sep = "\t", header = TRUE, fill = TRUE, quote = "", stringsAsFactors = TRUE)

MetabolicPaths_modules_select <- MetabolicPaths_modules %>% 
  select(bin_name, module,pathwise_module_completeness)


# Step 1: Create a complete data frame with all combinations of genome_name and function.
complete_MetabolicPaths_modules_select <- MetabolicPaths_modules_select %>%
  tidyr::expand(
    bin_name,
    module,
  )

complete_MetabolicPaths_modules_select <- complete_MetabolicPaths_modules_select %>%
  left_join(MetabolicPaths_modules_select, by = c("bin_name", "module"))

############################

# load genome metadata 
habitat_prefs <- read.table("/Users/home/SPECIES_LEVEL_PANGENOMES/TESTS/Classifying_genomes_to_oral_sites/Groups/Without_HP/Refined_habitat_prefs.txt", sep = "\t", header = TRUE, fill = TRUE, quote = "", stringsAsFactors = TRUE)  # Make sure to set the correct separator if it's not tab.

habitat_prefs <- habitat_prefs %>% 
  mutate(Genus = case_when(grepl('N_', Group_ID) ~ 'Neisseria',
                           grepl('K_', Group_ID) ~ 'Kingella',
                           grepl('E_', Group_ID) ~ 'Eikenella',
                           grepl('S_', Group_ID) ~ 'Simonsiella')) %>% 
  mutate(Species = case_when(grepl('E_corrodens', Group_ID) ~ 'corrodens',
                             grepl('E_exigua', Group_ID) ~ 'exigua',
                             grepl('E_glucosivorans', Group_ID) ~ 'glucosivorans',
                             grepl('E_halliae', Group_ID) ~ 'halliae',
                             grepl('K_bonacorsii_oralis', Group_ID) ~ 'K_oralis',
                             grepl('K_denitrificans', Group_ID) ~ 'denitrificans',
                             grepl('K_kingae', Group_ID) ~ 'kingae',
                             grepl('K_negevensis', Group_ID) ~ 'negevensis',
                             grepl('K_potus', Group_ID) ~ 'potus',
                             grepl('N_maigaei', Group_ID) ~ 'maigaei',
                             grepl('K_sp_str_SNUBH_2017', Group_ID) ~ 'sp_str_SNUBH_2017',
                             grepl('N_bacilliformis', Group_ID) ~ 'bacilliformis',
                             grepl('N_basseii', Group_ID) ~ 'basseii',
                             grepl('N_benedictiae', Group_ID) ~ 'benedictiae',
                             grepl('N_bergeri', Group_ID) ~ 'bergeri',
                             grepl('N_blantyrii', Group_ID) ~ 'blantyrii',
                             grepl('N_brasiliensis', Group_ID) ~ 'brasiliensis',
                             grepl('N_cinerea', Group_ID) ~ 'cinerea',
                             grepl('N_cinerea_str_CCUG_5746', Group_ID) ~ 'cinerea_str_CCUG_5746',
                             grepl('N_dentiae', Group_ID) ~ 'dentiae',
                             grepl('N_dumasiana', Group_ID) ~ 'dumasiana',
                             grepl('N_elongata', Group_ID) ~ 'elongata',
                             grepl('N_gonorrhoeae', Group_ID) ~ 'gonorrhoeae',
                             grepl('N_lactamica', Group_ID) ~ 'lactamica',
                             grepl('N_lactamica_str_NS19', Group_ID) ~ 'lactamica_str_NS19',
                             grepl('N_mucosa_', Group_ID) ~ 'mucosa',
                             grepl('N_oralis', Group_ID) ~ 'N_oralis',
                             grepl('N_polysaccharea', Group_ID) ~ 'polysaccharea',
                             grepl('N_shayeganii', Group_ID) ~ 'shayeganii',
                             grepl('N_sp_HMT_020', Group_ID) ~ 'sp_HMT_020',
                             grepl('N_sp_str_3986_and_str_51_81', Group_ID) ~ 'sp_str_3986_and_str_51_81',
                             grepl('N_sp_str_HSC_16F19', Group_ID) ~ 'sp_str_HSC_16F19',
                             grepl('N_sp_str_MVDL20_010259', Group_ID) ~ 'sp_str_MVDL20_010259',
                             grepl('N_subflava_I', Group_ID) ~ 'subflava_A',
                             grepl('N_subflava_II', Group_ID) ~ 'subflava_A',
                             grepl('N_subflava_III', Group_ID) ~ 'subflava_A',
                             grepl('N_subflava_IV', Group_ID) ~ 'subflava_A',
                             grepl('N_subflava_VII', Group_ID) ~ 'subflava_A',
                             grepl('N_subflava_VIII', Group_ID) ~ 'subflava_A',
                             grepl('N_subflava_V', Group_ID) ~ 'subflava_B',
                             grepl('N_uirgultaei', Group_ID) ~ 'uirgultaei',
                             grepl('N_viridiae', Group_ID) ~ 'viridiae',
                             grepl('N_wadsworthii', Group_ID) ~ 'wadsworthii',
                             grepl('N_weaveri', Group_ID) ~ 'weaveri',
                             grepl('N_meningitidis', Group_ID) ~ 'meningitidis',
                             grepl('S_muelleri', Group_ID) ~ 'muelleri'))




##### Remove undetected genomes from main df, which is any pangenome ID in the habitat_prefs df that has a granular preference of "NONE"
#make list of detected genomes
detected_genomes <- habitat_prefs %>% 
  filter(!(Habitat_Preference_Granular == "NONE")) %>% 
  select(item) %>% 
  droplevels()

genome_metadata <- habitat_prefs %>% 
  mutate(Pangenome_ID = item)

##############################


merged_df <- merge(complete_MetabolicPaths_modules_select, genome_metadata, by.x = "bin_name", by.y = "Pangenome_ID", all.x = TRUE)


# filter genomes to exclude undetected
merged_df <- merged_df %>% 
  filter(bin_name %in% detected_genomes$item) %>% 
  droplevels()


Modules_occurrence <- merged_df %>% 
  mutate(Occurrence = ifelse(pathwise_module_completeness >=0.75, 1, 0))

Modules_occurrence$Occurrence[is.na(Modules_occurrence$Occurrence)] <- 0


####  Make matrix for clustering ####

MetabolicPaths_modules_clustering_matrix <- Modules_occurrence %>% 
  select(bin_name, module, Occurrence)

# wide format
wide_df <- MetabolicPaths_modules_clustering_matrix %>% 
  pivot_wider(names_from = module,
              values_from = Occurrence)

wide_df <- as.data.frame(wide_df)

# rename row names to species IDs
wide_df2 <- wide_df[,-1]
rownames(wide_df2) <- wide_df[,1]

# convert to matrix
pathwise_module_completeness_matrix <- as.matrix(wide_df2)

# Convert NAs to zeros
pathwise_module_completeness_matrix[is.na(pathwise_module_completeness_matrix)] <- 0

# Drop columns with all zeros
pathwise_module_completeness_matrix <- pathwise_module_completeness_matrix[, colSums(pathwise_module_completeness_matrix) > 0]



# Create a mapping for multiple annotations
pathwise_module_completeness_annotation_row <- habitat_prefs %>%
  select(item, Habitat_Preference_Granular, Species, Genus) %>%
  rename(Genome = item, 'Habitat Preference' = Habitat_Preference_Granular) %>%
  filter(Genome %in% rownames(pathwise_module_completeness_matrix)) %>%
  column_to_rownames(var = "Genome")  # Set genome IDs as rownames

# Specify colors
pathwise_module_completeness_ann_colors = list(
  'Module Habitat Preference' = c(Plaque = "blue",Tongue = "red", 'Plaque and Gums' = "#25f4b9", None = "gray"),
  'Habitat Preference' = c(PLAQUE = "blue", TONGUE= "red", GUMS = "green"),
  Species = c(bacilliformis = "#c41c00", cinerea = "#06ea16", corrodens = "#01278f",denitrificans = "#0020c2",
              elongata = "#70eafa",exigua = "#f7e702",halliae = "#c2045d",mucosa = "#914747",
              muelleri = "#9a5bb5",N_oralis = "#ff00dd",sp_HMT_020 = "#8200f4",sp_str_SNUBH_2017 = "#f77cd2",
              subflava_A = "#ff006f", subflava_B = "orange", K_oralis = "red"),
  Genus = c(Eikenella = "#6363ea", Kingella = "#25f4b9", Neisseria = "#f46c6c", Simonsiella = "#edcb3b")
)



binary_colors <- c("black", "cyan")



# Create a lookup list for module habitats
module_habitat_lookup <- list(
  "Plaque" = c("M00880", "M00899", "M00554", "M00632", "M00012"),
  "Tongue" = c("M00026", "M00028", "M00015", "M00855"),
  "Plaque and Gums" = c("M00529")
)

# Extract all the unique modules from the Modules_occurrence file
all_modules <- unique(as.character(Modules_occurrence$module))

# Assign habitats to modules
module_habitat_df <- data.frame(
  module = all_modules,
  Habitat = "None", # Default habitat is "None"
  stringsAsFactors = FALSE
)

# Loop over each habitat and its corresponding modules
for (habitat in names(module_habitat_lookup)) {
  module_habitat_df$Habitat[module_habitat_df$module %in% module_habitat_lookup[[habitat]]] <- habitat
}

print(module_habitat_df)

# Create a mapping for multiple annotations
pathwise_module_completeness_annotation_col <- module_habitat_df %>%
  rename('Module Habitat Preference' = Habitat) %>% 
  column_to_rownames(var = "module") 



# Step 1: Create initial pheatmap to extract clustering information
pathwise_module_completeness_initial_plot <- pheatmap(
  pathwise_module_completeness_matrix,
  cluster_rows = TRUE,      
  cluster_cols = TRUE, 
  clustering_distance_rows = "binary",
  clustering_distance_cols = "binary",
  clustering_method = "complete",
  color = binary_colors, 
  show_rownames = FALSE,     # Hide row names if desired
  show_colnames = FALSE,
  annotation_row = pathwise_module_completeness_annotation_row, # Hide column names if desired
  annotation_colors = pathwise_module_completeness_ann_colors,
  annotation_col = pathwise_module_completeness_annotation_col,
  legend = FALSE,
  border_color = "black",
  annotation_legend = FALSE,
  treeheight_col = 50
)


save_pheatmap_pdf(pathwise_module_completeness_initial_plot, paste0(MainDIR,"20_FUNCTIONAL_ENRICHMENT/V2_KEGG_module_binary_pathway_completeness_75.pdf"), width = 9, height = 7)





```

```{r}


KEGG_occurence_habitat_binary_long <- KEGG_occurence_habitat_binary %>% 
  pivot_longer(
    cols = c(-KEGG_function, -KEGG_function_seq),
    names_to = "Habitat",
    values_to = "Presence"
  )



KEGG_occurence_long_summary <- KEGG_occurence_long %>%
  group_by(Habitat) %>% 
  summarize(n = sum(Presence_binary)/246)


KEGG_occurence_habitat_binary_long %>% 
  group_by(Habitat) %>% 
  summarize(Complete_Pathways = sum(Presence)) 




```


# 18. Gene-level analyses

Now that we have a list of potentially important gene clusters that are associated with each habitat, we can visualize their coverage.

Get gene cluster and gene ID info from pan DB.
```{bash, eval = FALSE}

# Get gene cluster presence/absence file from Pangenome DB
PAN="/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/09_PANGENOME/G_0213/P_0003_Neisseriaceae-RESULTS/P_0003_Neisseriaceae-PAN.db"
GENOMES="/Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/09_PANGENOME/G_0213/P_0003_Neisseriaceae-GENOMES.db"

sqlite3 $PAN ".mode tabs" ".headers on" "SELECT * FROM gene_clusters;" > /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/09_PANGENOME/G_0213/gene_clusters.tsv

sqlite3 $GENOMES ".mode tabs" ".headers on" "SELECT * FROM gene_info;" > /Users/home/SPECIES_LEVEL_PANGENOMES/Neisseriaceae/09_PANGENOME/G_0213/gene_info.tsv

```

The gene_info.tsv file contains: genome_name, gene_caller_id, aa_sequence, dna_sequence, partial and length.

The gene_clusters.tsv file contains: gene_caller_id, gene_cluster_id, genome_name and alignment_summary. 

The goal is to combine these two data frames with the various functional enrichment results above so that the functional enrichment results also have gene caller ID info and both aa_sequence and dna_sequence.


### Gene-level detection plots

*See Neisseriacea_gene_level_figures.RMD*


