#!/bin/bash


# variables
MAPQ=$1
projectID=P_0622_Haemophilus_Aggregatibacter
mainDIR=/workspace/jmarkwelchlab/$projectID
DIR_ContigsDB=04_CONTIGS_DB
DIR_Mapping=05_MAPPING/TEST_RUN/MAPQ_TESTS/$MAPQ
DIR_SinglePROF=06_SINGLE_PROFILE/TEST_RUN/MAPQ_TESTS/$MAPQ
DIR_MergedPROF=07_MERGED_PROFILE/TEST_RUN/MAPQ_TESTS/$MAPQ
minContigSIZE=200
threadsRun=2

CONTIGS_DB=$DIR_ContigsDB/${projectID}-contigs.db

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


# make merged profile for 3 oral sites
num_processes=3

declare -a siteArray=("BM" "TD" "PP")

for site in ${siteArray[@]}
do
	if [ ! -f "$DIR_MergedPROF/${MAPQ}_${site}/PROFILE.db" ]
	then
	((i=i%num_processes)); ((i++==0)) && wait
	/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/SCRIPTS/script-merged_profile.sh ${MAPQ} ${site} &
	fi
done

wait

# END merged profiling time
ELAPSED="${MAPQ}-Merged_profile $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> ${mainDIR}/time-00-overall.tsv
