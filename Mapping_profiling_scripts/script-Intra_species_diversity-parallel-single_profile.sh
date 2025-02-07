#!/bin/bash

# variables
projectID=$1
genome=$2
G_ID=$3
mainDIR=/workspace/jmarkwelchlab/$projectID
threadsRun=5
sampleList=$mainDIR/27_VARIABILITY/Intra_species_diversity/${genome}_samples.txt
DIR_Mapping=$mainDIR/27_VARIABILITY/Intra_species_diversity/05_MAPPING/${genome}
DIR_ContigsDB=$mainDIR/27_VARIABILITY/Intra_species_diversity/04_CONTIGS_DB
DIR_SinglePROF=$mainDIR/27_VARIABILITY/Intra_species_diversity/06_SINGLE_PROFILE/${genome}
CONTIGS_DB=$DIR_ContigsDB/${genome}-contigs.db

minContigSIZE=200

# change to directory
cd $mainDIR

# export variables
export mainDIR=$mainDIR
export threadsRun=$threadsRun
export minContigSIZE=$minContigSIZE
export DIR_Mapping=$DIR_Mapping
export DIR_SinglePROF=$DIR_SinglePROF

export CONTIGS_DB=$CONTIGS_DB

# START single profiling time
SECONDS=0

num_processes=6

for metagenome in `cat $sampleList`
do
        if [ ! -f "$DIR_SinglePROF/$metagenome/PROFILE.db" ]
        then
        ((i=i%num_processes)); ((i++==0)) && wait
        /workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/SCRIPTS/script-Intra_species_diversity-single_profile.sh $metagenome &
        fi
done

wait

# END single profiling time
ELAPSED="Single_profile ${genome} $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> ${mainDIR}/27_VARIABILITY/Intra_species_diversity/time-00-overall.tsv
