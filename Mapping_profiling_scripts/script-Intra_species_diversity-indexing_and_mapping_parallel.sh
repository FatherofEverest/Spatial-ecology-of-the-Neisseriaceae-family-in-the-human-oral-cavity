#!/bin/bash

# variables
projectID=$1
genome=$2
G_ID=$3
mainDIR=/workspace/jmarkwelchlab/$projectID
PROJECT_CONTIGS=$mainDIR/03_GENOMES_EDITED/${G_ID}.fa
threadsRun=10
sampleList=$mainDIR/27_VARIABILITY/Intra_species_diversity/${genome}_samples.txt
DIR_Mapping=$mainDIR/27_VARIABILITY/Intra_species_diversity/05_MAPPING/${genome}
DIR_Reads=$mainDIR/00_READS

# change to directory
cd $mainDIR

# generate REFERNCE-index for mapping
REFERENCE=$PROJECT_CONTIGS
INDEX=$DIR_Mapping/$G_ID

if [ ! -f "$INDEX.rev.2.bt2" ]
then
# START indexing time
SECONDS=0
bowtie2-build $REFERENCE $INDEX --threads $threadsRun
# END indexing time
ELAPSED="Reference_index $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> $mainDIR/27_VARIABILITY/Intra_species_diversity/time-00-overall.tsv
fi

# export variables
export mainDIR=$mainDIR
export threadsRun=$threadsRun
export DIR_Mapping=$DIR_Mapping
export DIR_Reads=$DIR_Reads
export INDEX=$INDEX

# START mapping time
SECONDS=0

num_processes=6

for metagenome in `cat $sampleList`
do
        if [ ! -f "$DIR_Mapping/$metagenome.bai" ]
        then
        ((i=i%num_processes)); ((i++==0)) && wait
        /workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/SCRIPTS/script-Intra_species_diversity-mapping_parallel.sh $metagenome &
        fi
done

wait

# END mapping time
ELAPSED="Mapping ${genome} $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> $mainDIR/27_VARIABILITY/Intra_species_diversity/time-00-overall.tsv
