#!/bin/bash

# variables
projectID=$1
mainDIR=/workspace/jmarkwelchlab/$projectID
threadsRun=15
sampleList=missing_samples_id-QC_IDs.txt
minContigSIZE=300
DIR_ContigsDB=04_CONTIGS_DB
DIR_Mapping=05_MAPPING
DIR_SinglePROF=06_SINGLE_PROFILE

CONTIGS_DB=$DIR_ContigsDB/${projectID}-contigs.db

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
	/workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/Mapping_profiling_scripts/script-single_profile.sh $metagenome &
	fi
done

wait

# END single profiling time
ELAPSED="013-Missing_Single_profile $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> ${mainDIR}/time-00-overall.tsv
