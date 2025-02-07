#!/bin/bash


# variables
projectID=$1

mainDIR=/workspace/jmarkwelchlab/$projectID
DIR_Contigs=03_GENOMES_EDITED
DIR_ContigsDB=04_CONTIGS_DB
DIR_MergedPROF=07_MERGED_PROFILE
DIR_SummaryPROF=08_PROFILE_SUMMARY
COLLECTION=Genomes
threadsRun=15

CONTIGS_DB=$DIR_ContigsDB/${projectID}-contigs.db
DECOMPOSE=$DIR_Contigs/${projectID}.decompose_modified.tsv

# change to directory
cd $mainDIR

# export variables
export mainDIR=$mainDIR
export DIR_MergedPROF=$DIR_MergedPROF
export DIR_Mapping=$DIR_Mapping
export DIR_SummaryPROF=$DIR_SummaryPROF
export minContigSIZE=$minContigSIZE
export threadsRun=$threadsRun
export COLLECTION=$COLLECTION
export DECOMPOSE=$DECOMPOSE
export CONTIGS_DB=$CONTIGS_DB

# make merged profile for 9 oral sites
num_processes=9

declare -a siteArray=("BM" "TD" "PP" "PT" "TH" "KG" "PB" "SA" "HP")

for site in ${siteArray[@]}
do
	if [ ! -f "$DIR_SummaryPROF/${projectID}_${site}-profile/bins_summary.txt" ]
	then
	((i=i%num_processes)); ((i++==0)) && wait
	/workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/Mapping_profiling_scripts/script-decompose_summarize_profile.sh ${projectID} ${site} &
	fi
done

wait

# END decompose and summarize time
ELAPSED="015-Decompose_estimate_SCG_taxonomy_summarize $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> ${mainDIR}/time-00-overall.tsv
