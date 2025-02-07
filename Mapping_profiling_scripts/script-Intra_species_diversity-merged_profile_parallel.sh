#!/bin/bash

# variables
projectID=$1
genome=$2
mainDIR=/workspace/jmarkwelchlab/$projectID
threadsRun=5
sampleList=$mainDIR/27_VARIABILITY/Intra_species_diversity/${genome}_samples.txt
DIR_Mapping=$mainDIR/27_VARIABILITY/Intra_species_diversity/05_MAPPING/${genome}
DIR_ContigsDB=$mainDIR/27_VARIABILITY/Intra_species_diversity/04_CONTIGS_DB
CONTIGS_DB=$DIR_ContigsDB/${genome}-contigs.db
DIR_SinglePROF=$mainDIR/27_VARIABILITY/Intra_species_diversity/06_SINGLE_PROFILE/${genome}
DIR_MergedPROF=$mainDIR/27_VARIABILITY/Intra_species_diversity/07_MERGED_PROFILE/${genome}
minContigSIZE=200
threadsRun=5

# change to directory
cd $mainDIR

# export variables
export mainDIR=$mainDIR
export DIR_SinglePROF=$DIR_SinglePROF
export DIR_MergedPROF=$DIR_MergedPROF
export DIR_Mapping=$DIR_Mapping
export minContigSIZE=$minContigSIZE
export threadsRun=$threadsRun
export CONTIGS_DB=$CONTIGS_DB

# START merged profiling time
SECONDS=0


# make merged profile for 8 oral sites
num_processes=3

declare -a siteArray=("BM" "TD" "PP")

for site in ${siteArray[@]}
do
        if [ ! -f "$DIR_MergedPROF/${genome}_${site}/PROFILE.db" ]  
        then
        ((i=i%num_processes)); ((i++==0)) && wait
        /workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/SCRIPTS/script-Intra_species_diversity-merged_profile.sh ${genome} ${site} &
        fi
done

wait

# END merged profiling time
ELAPSED="Merged_profile ${genome} $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> $mainDIR/27_VARIABILITY/Intra_species_diversity/time-00-overall.tsv