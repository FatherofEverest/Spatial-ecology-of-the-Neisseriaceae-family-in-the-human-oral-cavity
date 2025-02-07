#!/bin/bash

# variables
projectID=$1
genome=$2
mainDIR=/workspace/jmarkwelchlab/$projectID

DIR_ContigsDB=$mainDIR/27_VARIABILITY/Intra_species_diversity/04_CONTIGS_DB
DIR_MergedPROF=$mainDIR/27_VARIABILITY/Intra_species_diversity/07_MERGED_PROFILE
DIR_Variability=$mainDIR/27_VARIABILITY/Intra_species_diversity/08_VARIABILITY

threadsRun=5


# change to directory
cd $mainDIR

# export variables
export mainDIR=$mainDIR
export DIR_MergedPROF=$DIR_MergedPROF
export DIR_ContigsDB=$DIR_ContigsDB
export DIR_Variability=$DIR_Variability
export threadsRun=$threadsRun

# START timer
SECONDS=0


# make variability reports for 3 oral sites
num_processes=3

declare -a siteArray=("BM" "TD" "PP")

for site in ${siteArray[@]}
do
        if [ ! -f "$DIR_Variability/${genome}/${site}/${genome}_SNVs.txt" ]  
        then
        ((i=i%num_processes)); ((i++==0)) && wait
        /workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/SCRIPTS/script-Intra_species_diversity-variability.sh ${genome} ${site} &
        fi
done

wait

# END merged profiling time
ELAPSED="Variability ${genome} $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> $mainDIR/27_VARIABILITY/Intra_species_diversity/time-00-overall.tsv









