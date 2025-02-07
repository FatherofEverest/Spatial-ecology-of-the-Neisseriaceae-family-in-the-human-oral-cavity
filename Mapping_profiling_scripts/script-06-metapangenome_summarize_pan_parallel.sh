#!/bin/bash


# variables
projectID=$2

# pangenome collection
panCOL=$1

mainDIR=/workspace/jmarkwelchlab/$projectID
DIR_ContigsDB=04_CONTIGS_DB
DIR_MergedPROF=07_MERGED_PROFILE
DIR_Pangenome=09_PANGENOME
DIR_SummaryPAN=10_PANGENOME_SUMMARY

PAN_LAYERS=$DIR_Pangenome/${projectID}-add_info.layers.tsv

# change to directory
cd $mainDIR

# export variables
export mainDIR=$mainDIR
export panCOL=$panCOL
export DIR_ContigsDB=$DIR_ContigsDB
export DIR_MergedPROF=$DIR_MergedPROF
export DIR_Pangenome=$DIR_Pangenome
export DIR_SummaryPAN=$DIR_SummaryPAN
export PAN_LAYERS=$PAN_LAYERS

# START metapangenome time
SECONDS=0

# make merged profile for 9 oral sites
num_processes=9

declare -a siteArray=("BM" "TD" "PP" "PT" "TH" "KG" "PB" "SA" "HP")

for site in ${siteArray[@]}
do
	if [ ! -f "$DIR_Pangenome/HMP_${site}/${projectID}-RESULTS/${projectID}-PAN.db" ]
	then
	((i=i%num_processes)); ((i++==0)) && wait
	/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/SCRIPTS/script-metapangenome.sh ${projectID} ${site} &
	fi
done

wait

# END metapangenome time
ELAPSED="017-Metapangenome $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> ${mainDIR}/time-00-overall.tsv

# START summary metapangenome time
SECONDS=0

# make merged profile for 9 oral sites
num_processes=10

declare -a siteArray=("BM" "TD" "PP" "PT" "TH" "KG" "PB" "SA" "HP")

for site in ${siteArray[@]}
do
	if [ ! -d "$DIR_SummaryPAN/HMP_${site}" ]
	then
	((i=i%num_processes)); ((i++==0)) && wait
	/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/SCRIPTS/script-summary_metapan.sh ${projectID} ${site} &
	fi
done

wait

# END summary metapangenome time
ELAPSED="018-Summary_metapangenome $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> ${mainDIR}/time-00-overall.tsv
