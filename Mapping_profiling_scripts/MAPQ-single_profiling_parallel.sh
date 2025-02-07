#!/bin/bash

# variables
MAPQ=$1
projectID=P_0622_Haemophilus_Aggregatibacter
mainDIR=/workspace/jmarkwelchlab/$projectID
threadsRun=2
sampleList=DATA/MAPQ_TEST_metagenome_IDs.txt
minContigSIZE=200
DIR_ContigsDB=04_CONTIGS_DB
DIR_Mapping=05_MAPPING/TEST_RUN/MAPQ_TESTS/${MAPQ}_MAPQ
DIR_SinglePROF=06_SINGLE_PROFILE/TEST_RUN/MAPQ_TESTS/${MAPQ}_MAPQ

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
	/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/SCRIPTS/script-single_profile.sh ${metagenome}_MAPQ_${i} &
	fi
done

wait

# END single profiling time
ELAPSED="MAPQ_TEST-Single_profile $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> ${mainDIR}/time-00-overall.tsv
