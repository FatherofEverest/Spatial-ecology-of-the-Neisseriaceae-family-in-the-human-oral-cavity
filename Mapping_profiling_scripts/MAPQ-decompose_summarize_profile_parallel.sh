#!/bin/bash


# variables
MAPQ=$1
projectID=P_0622_Haemophilus_Aggregatibacter
mainDIR=/workspace/jmarkwelchlab/$projectID
DIR_Contigs=03_GENOMES_EDITED
DIR_ContigsDB=04_CONTIGS_DB
DIR_MergedPROF=07_MERGED_PROFILE/TEST_RUN/MAPQ_TESTS/$MAPQ
DIR_SummaryPROF=08_PROFILE_SUMMARY/TEST_RUN/MAPQ_TESTS/$MAPQ
COLLECTION=Genomes
threadsRun=2

CONTIGS_DB=$DIR_ContigsDB/${projectID}-contigs.db
DECOMPOSE=$DIR_Contigs/${MAPQ}.decompose.tsv

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
num_processes=3

declare -a siteArray=("BM" "TD" "PP")

for site in ${siteArray[@]}
do
	if [ ! -f "$DIR_SummaryPROF/${MAPQ}_${site}-profile/bins_summary.txt" ]
	then
	((i=i%num_processes)); ((i++==0)) && wait
	/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/SCRIPTS/script-decompose_summarize_profile.sh ${MAPQ} ${site} &
	fi
done

wait

# END decompose and summarize time
ELAPSED="${MAPQ}-Decompose_estimate_SCG_taxonomy_summarize $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> ${mainDIR}/time-00-overall.tsv
