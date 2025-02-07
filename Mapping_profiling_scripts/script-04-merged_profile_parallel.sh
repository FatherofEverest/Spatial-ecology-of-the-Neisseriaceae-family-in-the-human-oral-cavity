#!/bin/bash


# variables
projectID=$1
mainDIR=/workspace/jmarkwelchlab/$projectID
DIR_ContigsDB=04_CONTIGS_DB
DIR_Mapping=05_MAPPING
DIR_SinglePROF=06_SINGLE_PROFILE
DIR_MergedPROF=07_MERGED_PROFILE
minContigSIZE=300
threadsRun=15

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

# make merged profile for HP (unique-sample)
site=HP

if [ ! -f "$DIR_MergedPROF/${projectID}_${site}/PROFILE.db" ]
then
/workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/Mapping_profiling_scripts/script-merged_profile-HP.sh ${projectID} ${site} &
fi

# make merged profile for 8 oral sites
num_processes=8

declare -a siteArray=("BM" "TD" "PP" "PT" "TH" "KG" "PB" "SA")

for site in ${siteArray[@]}
do
	if [ ! -f "$DIR_MergedPROF/${projectID}_${site}/PROFILE.db" ]
	then
	((i=i%num_processes)); ((i++==0)) && wait
	/workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/Mapping_profiling_scripts/script-merged_profile.sh ${projectID} ${site} &
	fi
done

wait

# END merged profiling time
ELAPSED="014-Merged_profile $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> ${mainDIR}/time-00-overall.tsv
