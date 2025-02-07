####################################################################################################################
##### download_genomes.sh
####################################################################################################################

#!/bin/bash

projectID=P_0003_Neisseriaceae
mainDIR=/workspace/jmarkwelchlab
DIR_NCBI=$mainDIR/$projectID/01_NCBI_GENOMES
downloadNCBI=$mainDIR/$projectID/DATA/03_P_0003_Neisseriaceae-download_url.txt

for gca_assembly in `cat $downloadNCBI`
do
gcaFile=$( echo $gca_assembly | awk -F'/' '{print $NF}' | sed -e 's/\.gz//' )
if [ ! -f $DIR_NCBI/${gcaFile} ]
then
echo "File not found - Downloading: ${gcaFile} "
wget -q -P $DIR_NCBI $gca_assembly
gunzip $DIR_NCBI/${gcaFile}.gz
echo "unzipped: ${gcaFile} "
fi
done

####################################################################################################################
##### rename_genomes.sh
####################################################################################################################

#!/bin/bash

projectID=P_0003_Neisseriaceae
mainDIR=/workspace/jmarkwelchlab
DIR_NCBI=$mainDIR/$projectID/01_NCBI_GENOMES
DIR_Assemblies=$mainDIR/$projectID/02_ASSEMBLIES
DIR_Data=$mainDIR/$projectID/DATA
projectMetadata=$DIR_Data/02_${projectID}.csv
genomesNCBIid=$DIR_Data/08_P_0003_Neisseriaceae-ncbi_genomes_id.txt
rawNCBIid=$DIR_Data/04_P_0003_Neisseriaceae-ncbi_raw_genomes_id.txt
genomesID=$DIR_Data/id_genomes.txt

# create a list of unique IDs containing a number and the assembly ID separated by "_id_"; use GCA assembly ID column 9
# (e.g. G_0001_id_GCA_000000000.0)
cat $projectMetadata | awk -F',' -v OFS="\t" 'NR>1{ printf "G_""%04i_id_%s\n", NR-1,$12 }' | sed -e 's/"//g'  > $rawNCBIid

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

####################################################################################################################
##### reformat_genome_deflines_and_build_contig_databases.sh
####################################################################################################################


#!/bin/bash
projectID=P_0003_Neisseriaceae
mainDIR=/workspace/jmarkwelchlab
DIR_Contigs=$mainDIR/$projectID/03_GENOMES_EDITED
DIR_ContigsDB=$mainDIR/$projectID/04_CONTIGS_DB
DIR_Assemblies=$mainDIR/$projectID/02_ASSEMBLIES
DIR_Data=$mainDIR/$projectID/DATA
genomesID=$DIR_Data/id_genomes.txt


for genome in `cat $genomesID`
do
EDITED_GENOMES=$DIR_Assemblies/${genome}-RAW.fa
CONTIGS=$DIR_Contigs/${genome}.fa
REPORT=$DIR_Contigs/${genome}.report.tsv
anvi-script-reformat-fasta -o $CONTIGS --simplify-names --prefix ${genome} --seq-type NT -r $REPORT $EDITED_GENOMES 
done



for genome in `cat $genomesID`
do
CONTIGS=$DIR_Contigs/${genome}.fa
contigsDB=$DIR_ContigsDB/${genome}-contigs.db
numThreads=5
anvi-gen-contigs-database -f $CONTIGS -n ${genome} -o $contigsDB -T $numThreads 
done


####################################################################################################################
##### anvio_hmms_trnas_gtdb.sh
####################################################################################################################


#!/bin/bash

projectID=P_0003_Neisseriaceae
mainDIR=/workspace/jmarkwelchlab
DIR_ContigsDB=$mainDIR/$projectID/04_CONTIGS_DB
DIR_Data=$mainDIR/$projectID/DATA
genomesID=$DIR_Data/id_genomes.txt

N=20
(
for genome in `cat $genomesID`
do
((i=i%N)); ((i++==0)) && wait
contigsDB=$DIR_ContigsDB/${genome}-contigs.db
numThreads=5
echo "Running anvio hmms for ${genome}"
anvi-run-hmms -c $contigsDB -T $numThreads &
done
)


# search for tRNAs in genomes (scan-trnas)
N=20
(
for genome in `cat $genomesID`
do
((i=i%N)); ((i++==0)) && wait
contigsDB=$DIR_ContigsDB/${genome}-contigs.db
numThreads=5
echo "Running anvio scan trnas for ${genome}"
anvi-scan-trnas -c $contigsDB -T $numThreads &
done
)

# GTDB taxonomy 

module load diamond

for genome in `cat $genomesID`
do
contigsDB=$DIR_ContigsDB/${genome}-contigs.db
numThreads=5
echo "Running anvio scg taxonomy for ${genome}"
anvi-run-scg-taxonomy -c $contigsDB -T $numThreads --write-buffer-size 30000 

####################################################################################################################
##### check_anvio_contigDB_hmms.sh
####################################################################################################################


#!/bin/bash
projectID=P_0003_Neisseriaceae
mainDIR=/workspace/jmarkwelchlab
DIR_ContigsDB=$mainDIR/$projectID/04_CONTIGS_DB
contigsDB=$DIR_ContigsDB

for contig in `ls $contigsDB`
  do
      anvi-display-contigs-stats  --report-as-text -o $mainDIR/$projectID/19_Contig_db_stats/$contig-stats  $DIR_ContigsDB/$contig
done

# check each contig-stats file 
echo -e "contig_db\tBacteria_71" > $mainDIR/$projectID/19_Contig_db_stats/Contig_bacteria_71_summary.txt

for contig in `ls $mainDIR/$projectID/19_Contig_db_stats`
  do
    echo "##### STARTING WITH $contig #####"
    cat 19_Contig_db_stats/$contig | grep 'contigs_db\|Bacteria_71' | sed 'n; n; d' > tmp
    contig_db=$(awk 'NR==1{print $2}' tmp)
    Bacteria_71=$(awk 'NR==2{print $2}' tmp)
    echo -e "$contig_db\t$Bacteria_71" >> $mainDIR/$projectID/19_Contig_db_stats/Contig_bacteria_71_summary.txt
    rm tmp
    echo "##### DONE WITH $contig #####"
done


####################################################################################################################
##### Rename_HOMD_create_ITEMS_Contigs_path_file.sh
####################################################################################################################


#!/bin/bash

projectID=P_0003_Neisseriaceae
mainDIR=/workspace/jmarkwelchlab
DIR_Data=$mainDIR/$projectID/DATA
nameConversions=$DIR_Data/13_${projectID}-name_conversions.txt
rawNCBIid=$DIR_Data/new_04_${projectID}-ncbi_raw_genomes_id.txt
seqidHOMD=/workspace/jmarkwelchlab/HOMD_INFO/SEQID_info.txt
newHOMDid=$DIR_Data/07_${projectID}-homd_id.txt
ITEMS=${projectID}-add_info.items.txt
projectMetadata=$DIR_Data/02_${projectID}.csv
genomesID=$DIR_Data/id_genomes.txt

# get conversion names
# head is assembly_accession,host,isolation_source,Genus,Species,Strain_abrv,ftp_path
while IFS= read -r line
do
genomeNo=$( echo $line | awk -F'_id_' '{print $1}')
assemblyID=$( echo $line | awk -F'_id_' '{print $2}')
putGenoName=$( grep "$assemblyID" $projectMetadata | awk -F',' '{gsub(/ /, "_", $6); printf "%s_%s_str_%s_id_%s\n", substr($4,1,1), $5, $6, $1}' | sed -e 's/\./_/')
echo -e "$putGenoName\t$genomeNo" >> $nameConversions
done < $rawNCBIid

# find HOMD ID
for hypo_hmt in `cat $rawNCBIid | awk -F'_id_' '{print $2}' `
do
cat $seqidHOMD | grep "$hypo_hmt" | awk -v gca_id="$hypo_hmt" 'BEGIN{FS=OFS="\t"}{print $3"_"$4"_str_"$5"_id_"gca_id}' | sed -e's/ //g' | sed -e 's/_sp\./_sp_/' | awk 'BEGIN{FS=OFS="_"}{ $1= substr($1,1,1)}1' | awk '{if ($1 ~/HMT/) gsub("_sp","");print}' | awk -F'_id_' -v OFS="\t" '{print $NF,$0}' | sed -e 's/-/_/g' | grep 'HMT' | sed -e 's/HMT/sp_HMT_/' >> $newHOMDid
done

# edit name conversion with HOMD sps.
while IFS= read -r line
do
gcaID=$( echo $line | awk '{print $1}' | sed -e 's/\./_/' )
oldID=$( grep "$gcaID" $nameConversions | awk '{print $1}' )
newID=$( echo $line | awk '{print $2}' | sed -e 's/\./_/' )
sed -i "s/$oldID/$newID/" $nameConversions
done < $newHOMDid


# make meta data file for genome renaming
echo -e "item\tGenome_ID\tSpecies\tAssembly_ID\tAssembly_ID_norm\tStrain\tGenome_in_HOMD\tHOMD_ID\tG_ID" > $ITEMS
for gcaID in `cat $nameConversions | cut -f1 | awk -F'_id_' '{print $2}' | rev | sed -e 's/_/./' | rev `
do
gcaID_ed=$(echo $gcaID | sed -e 's/\./_/')
itemsID=$(grep "$gcaID_ed" $nameConversions | cut -f1 | awk -F'_str_' -v OFS="\t" '{print $0,$0,$1}')
strainID=$( grep "$gcaID_ed" $nameConversions | cut -f1 | awk -F'_str_' -v OFS="\t" '{print $2}' | awk -F'_id_' -v OFS="\t" '{print $1}' )
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
echo -e "$itemsID\t$gcaID_ed\t$gcaID\t$strainID\t$inHOMD\t$hmtID\t$gID" >> $ITEMS
done


# Contig path file
GENOMES=$mainDIR/$projectID/DATA/$projectID-contig_paths.txt
echo -e 'name\tcontigs_db_path' > $GENOMES
for genome in `cat $genomesID`
do
cat $nameConversions | grep "$genome" | awk -v workPath="$PWD" 'BEGIN{FS=OFS="\t"}{print $1,workPath"/04_CONTIGS_DB/"$2"-contigs.db"}' >> $GENOMES
done



####################################################################################################################
##### decontaminated_pangenome.sh
####################################################################################################################


#!/bin/bash

export BLASTDB_LMDB_MAP_SIZE=1000000000

projectID=P_0003_Neisseriaceae
mainDIR=/workspace/jmarkwelchlab
DIR_Pangenome=$mainDIR/$projectID/09_PANGENOME
PanProject=${projectID}-decontaminated-pangenome
PanDir=$DIR_Pangenome/$PanProject
GENOMES_decontaminated=$mainDIR/$projectID/DATA/${projectID}-decontaminated-contig_paths.txt
GENOMES_decontaminated_pangenome_DB=$PanDir/${PanProject}-GENOMES.db
PAN_DIR=$PanDir/${PanProject}-RESULTS
PAN_DB=$PanDir/${PanProject}-RESULTS/${PanProject}-PAN.db
ANI_DIR=$PanDir/${PanProject}-RESULTS/ANI_RESULTS
layersADD_decontaminated=DATA/$projectID-decontaminated-add_info.items.txt
PAN_DES_decontaminated=DATA/description_decontaminated_pangenome.txt
TH=10


# external Genome storage
anvi-gen-genomes-storage -e $GENOMES_decontaminated -o $GENOMES_decontaminated_pangenome_DB

# external pangenome
mu -g $GENOMES_decontaminated_pangenome_DB --use-ncbi-blast --align-with muscle  --minbit 0.5 --mcl-inflation 10 -n ${PanProject} -o $PAN_DIR --num-threads $TH --enforce-hierarchical-clustering --description $PAN_DES_decontaminated

# add LAYERS
anvi-import-misc-data -t layers -p $PAN_DB $layersADD_decontaminated

# external genome similarity
anvi-compute-genome-similarity -e $GENOMES_decontaminated -o $ANI_DIR -p $PAN_DB --program pyANI --method ANIb --num-threads $TH --log-file initial_pangenome_pyANIlog




####################################################################################################################
##### initial_decontaminated_pangenome_single_copy_core_gene_phylogeny.sh
####################################################################################################################


#!/bin/bash


projectID=P_0003_Neisseriaceae
mainDIR=/workspace/jmarkwelchlab
DIR_Pangenome=$mainDIR/$projectID/09_PANGENOME
PanProject=${projectID}-decontaminated-pangenome
PanDir=$DIR_Pangenome/$PanProject
GENOMES_decontaminated_pangenome_DB=$PanDir/${PanProject}-GENOMES.db
PAN_DIR=$PanDir/${PanProject}-RESULTS
PAN_DB=$PanDir/${PanProject}-RESULTS/${PanProject}-PAN.db
N_genomes=473



anvi-get-sequences-for-gene-clusters -g $GENOMES_decontaminated_pangenome_DB \
                                     -p $PAN_DB \
                                     -o $PAN_DIR/single_copy_core_genes-fasta \
                                     --max-num-genes-from-each-genome 1 \
                                     --min-num-genomes-gene-cluster-occurs $N_genomes \
                                     --min-geometric-homogeneity-index 1 \
                                     --concatenate-gene-clusters 
                                     
                                     

iqtree -s $PAN_DIR/single_copy_core_genes-fasta -nt AUTO -m WAG -bb 1000                                      
                                     
                                    


####################################################################################################################
##### contig_lengths_per_genome.sh
####################################################################################################################



#!/bin/bash

projectID=P_0003_Neisseriaceae
mainDIR=/workspace/jmarkwelchlab
DIR_Assemblies=$mainDIR/$projectID/02_ASSEMBLIES
DIR_OUT=$mainDIR/$projectID/CONTIG_LENGTH_TEST
DIR_Data=$mainDIR/$projectID/DATA
genomesID=$DIR_Data/id_genomes.txt


echo -e "genome\tcontig_length_threshold\tcontigs_removed\tTotal_no_contigs\tpercent_contigs_removed" > CONTIG_LENGTH_TEST/final_table.txt

for length in $(seq 100 100 1000); do 
   for genome in `cat $genomesID`; do
    anvi-script-reformat-fasta -l $length \
      -o $DIR_OUT/${genome}_${length}.fa \
      --simplify-names \
      --prefix $genome \
      --seq-type NT \
      $DIR_Assemblies/$genome-RAW.fa 2>&1 | tee $DIR_OUT/${genome}_${length}_terminal_output.txt 
    
    contigs_removed=$(cat $DIR_OUT/${genome}_${length}_terminal_output.txt | grep "Contigs removed" | sed -e 's/^.*: //g' | sed -e 's/ .*//g') 
    Total_no_contigs=$(cat $DIR_OUT/${genome}_${length}_terminal_output.txt | grep "Total num contigs" | sed -e 's/^.* //g') 
    percent_contigs_removed=$(cat $DIR_OUT/${genome}_${length}_terminal_output.txt | grep "Contigs removed" | sed -e 's/^.*(//g' | sed -e 's/ of all)//g') 

    echo -e "$genome\t$length\t$contigs_removed\t$Total_no_contigs\t$percent_contigs_removed" >> $DIR_OUT/final_table.txt 
    
    rm $DIR_OUT/${genome}_${length}.fa $DIR_OUT/${genome}_${length}_terminal_output.txt 
   done 
done




####################################################################################################################
##### contig_stats_per_genome.sh
####################################################################################################################


#!/bin/bash

projectID=P_0003_Neisseriaceae
mainDIR=/workspace/jmarkwelchlab
DIR_Contigs=$mainDIR/$projectID/04_CONTIGS_DB
DIR_OUT=$mainDIR/$projectID/CONTIG_LENGTH_TEST
DIR_Data=$mainDIR/$projectID/DATA
genomesID=$DIR_Data/id_genomes.txt


echo -e "Genome\tLongest_Contig\tShortest_Contig\tNum_Contigs_Greater_100kb\tNum_Contigs_Greater_50kb\tNum_Contigs_Greater_20kb\tNum_Contigs_Greater_10kb\tNum_Contigs_Greater_5kb\tNum_Contigs_Greater_2_5kb" > $DIR_OUT/Contig_stats.txt


for genome in `cat $genomesID`; do
  anvi-display-contigs-stats $DIR_Contigs/${genome}-contigs.db --report-as-text  -o $DIR_OUT/${genome}-contigsDB_stats.txt
  
  Longest_Contig=$(cat $DIR_OUT/${genome}-contigsDB_stats.txt | grep "Longest Contig" | cut -f2 )
  Shortest_Contig=$(cat $DIR_OUT/${genome}-contigsDB_stats.txt | grep "Shortest Contig" | cut -f2 )
  Num_Contigs_Greater_100kb=$(cat $DIR_OUT/${genome}-contigsDB_stats.txt | grep "Num Contigs > 100 kb" | cut -f2 )
  Num_Contigs_Greater_50kb=$(cat $DIR_OUT/${genome}-contigsDB_stats.txt | grep "Num Contigs > 50 kb" | cut -f2 )
  Num_Contigs_Greater_25kb=$(cat $DIR_OUT/${genome}-contigsDB_stats.txt | grep "Num Contigs > 25 kb" | cut -f2 )
  Num_Contigs_Greater_10kb=$(cat $DIR_OUT/${genome}-contigsDB_stats.txt | grep "Num Contigs > 10 kb" | cut -f2 )
  Num_Contigs_Greater_5kb=$(cat $DIR_OUT/${genome}-contigsDB_stats.txt | grep "Num Contigs > 5 kb" | cut -f2 )
  Num_Contigs_Greater_2_5kb=$(cat $DIR_OUT/${genome}-contigsDB_stats.txt | grep "Num Contigs > 2.5 kb" | cut -f2 )
  
  echo -e "$genome\t$Longest_Contig\t$Shortest_Contig\t$Num_Contigs_Greater_100kb\t$Num_Contigs_Greater_50kb\t$Num_Contigs_Greater_25kb\t$Num_Contigs_Greater_10kb\t$Num_Contigs_Greater_5kb\t$Num_Contigs_Greater_2_5kb" >> $DIR_OUT/Contig_stats.txt
  
  rm $DIR_OUT/${genome}-contigsDB_stats.txt

done






#####################

#!/bin/bash

# variables
projectID=P_0003_Neisseriaceae
mainDIR=/workspace/jmarkwelchlab
DIR_CLE=$mainDIR/$projectID/CONTIG_LENGTH_TEST
DIR_Contigs=$DIR_CLE/03_GENOMES_EDITED
DIR_Mapping=$DIR_CLE/05_MAPPING
EDITED_CONCAT_GENOMES_BOTH=$DIR_Contigs/Edited_Concatenated_both.fa
DIR_Reads=$mainDIR/P_0622_Haemophilus_Aggregatibacter/00_READS
sampleList=$DIR_CLE/Top_30_TD_samples.txt
threadsRun=10

# change to directory
cd $mainDIR

# generate REFERNCE-index for mapping
REFERENCE=$EDITED_CONCAT_GENOMES_BOTH
INDEX=$DIR_Mapping/Concatenated_reference_set

if [ ! -f "$INDEX.rev.2.bt2" ]
then
# START indexing time
SECONDS=0
bowtie2-build $REFERENCE $INDEX --threads $threadsRun
# END indexing time
ELAPSED="006-Reference_index $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> ${mainDIR}/time-00-overall.tsv
fi

# export variables
export mainDIR=$mainDIR
export threadsRun=$threadsRun
export DIR_Mapping=$DIR_Mapping
export DIR_Reads=$DIR_Reads
export INDEX=$INDEX

# START mapping time
SECONDS=0

num_processes=10

for metagenome in `cat $sampleList | awk 'NR>1{print $4}'`
do
        if [ ! -f "$DIR_Mapping/$metagenome.bai" ]
        then
        ((i=i%num_processes)); ((i++==0)) && wait
        /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/Contig_length_test_script-mapping_parallel.sh $metagenome &
        fi
done

wait

# END mapping time
ELAPSED="Mapping $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> ${DIR_CLE}/time.tsv





##########  /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/Contig_length_test_script-mapping_parallel.sh


#!/bin/bash

sample=$1

# START mapping time per sample
SECONDS=0

# readsVAR
READS=$DIR_Reads/$sample

# mappingVAR
MAP=$DIR_Mapping/$sample

# mapping Bowtie2 (very-sensitive, end-to-end, no-unal)
bowtie2 -x $INDEX -1 ${READS}_R1.fastq.gz -2 ${READS}_R2.fastq.gz -S $MAP.sam --end-to-end --very-sensitive --no-unal --threads $threadsRun

# SAM to BAM
samtools view -F 4 -bS $MAP.sam -o ${MAP}-RAW.bam --threads $threadsRun

# sort BAM
samtools sort -o $MAP.bam ${MAP}-RAW.bam --threads $threadsRun

# BAM index
samtools index $MAP.bam $MAP.bai

# remove intermediate files
rm ${MAP}-RAW.bam $MAP.sam

# END mapping time per sample
ELAPSED="${sample} mapping time $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> ${DIR_CLE}/time.tsv





#!/bin/bash

# variables
projectID=P_0003_Neisseriaceae
mainDIR=/workspace/jmarkwelchlab/$projectID
DIR_CLE=$mainDIR/$projectID/CONTIG_LENGTH_TEST
DIR_ContigsDB=$DIR_CLE/04_CONTIGS_DB
DIR_Mapping=$DIR_CLE/05_MAPPING
DIR_SinglePROF=$DIR_CLE/06_PROFILES

CONTIGS_DB=$DIR_ContigsDB/Concatenated_both-contigs.db
sampleList=$DIR_CLE/Top_30_TD_samples.txt
threadsRun=10

# change to directory
cd $mainDIR

# export variables
export mainDIR=$mainDIR
export threadsRun=$threadsRun
export DIR_Mapping=$DIR_Mapping
export DIR_SinglePROF=$DIR_SinglePROF

export CONTIGS_DB=$CONTIGS_DB

# START single profiling time
SECONDS=0

num_processes=10

for metagenome in `cat $sampleList | awk 'NR>1{print $4}'`
do
        if [ ! -f "$DIR_SinglePROF/$metagenome/PROFILE.db" ]
        then
        ((i=i%num_processes)); ((i++==0)) && wait
        /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/contig_length_test_script-single_profile.sh $metagenome &
        fi
done

wait

# END single profiling time
ELAPSED="Profiling $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> ${DIR_CLE}/time.tsv



#!/bin/bash

sample=$1

# mappingVAR
BAM=$DIR_Mapping/$sample.bam

# single profle directory
SINGLE_PROFILE=$DIR_SinglePROF/$sample

# single profiling per sample
anvi-profile -i $BAM -c $CONTIGS_DB -o $SINGLE_PROFILE -S $sample --num-threads $threadsRun --write-buffer-size-per-thread 30000




