#!/bin/bash

# variables
projectID=P_0003_Neisseriaceae
mainDIR=/workspace/jmarkwelchlab
DIR_CLE=$mainDIR/$projectID/CONTIG_LENGTH_TEST
DIR_ContigsDB=$DIR_CLE/04_CONTIGS_DB
DIR_Mapping=$DIR_CLE/05_MAPPING
DIR_SinglePROF=$DIR_CLE/06_PROFILES/Concatenated_both_short_1000_filter

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
export DIR_CLE=$DIR_CLE
export CONTIGS_DB=$CONTIGS_DB

# START single profiling time
SECONDS=0

num_processes=10

for metagenome in `cat $sampleList | awk 'NR>1{print $4}'`
do
        if [ ! -f "$DIR_SinglePROF/$metagenome/PROFILE.db" ]
        then
        ((i=i%num_processes)); ((i++==0)) && wait
        /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/contig_length_test_script-single_profile_short_1000_filter.sh $metagenome &
        fi
done

wait

# END single profiling time
ELAPSED="Profiling - short 1000 nt filter $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> ${DIR_CLE}/time.tsv

