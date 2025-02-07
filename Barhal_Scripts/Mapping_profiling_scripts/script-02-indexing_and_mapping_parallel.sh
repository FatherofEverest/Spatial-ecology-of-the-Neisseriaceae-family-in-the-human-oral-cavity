#!/bin/bash

# variables
projectID=$1
mainDIR=/workspace/jmarkwelchlab/$projectID
PROJECT_CONTIGS=03_GENOMES_EDITED/${projectID}.fa
threadsRun=15
sampleList=samples_id-QC_IDs.txt
DIR_Mapping=05_MAPPING
DIR_Reads=00_READS

# change to directory
cd $mainDIR

# generate REFERNCE-index for mapping
REFERENCE=$PROJECT_CONTIGS
INDEX=$DIR_Mapping/$projectID

if [ ! -f "$INDEX.rev.2.bt2" ]
then
# START indexing time
SECONDS=0
bowtie2-build $REFERENCE $INDEX --threads $threadsRun
# END indexing time
ELAPSED="006-Reference_index $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> ${mainDIR}/time-00-overall.tsv
fi

# export variables
export mainDIR=$mainDIR
export threadsRun=$threadsRun
export DIR_Mapping=$DIR_Mapping
export DIR_Reads=$DIR_Reads
export INDEX=$INDEX

# START mapping time
SECONDS=0

num_processes=6

for metagenome in `cat $sampleList`
do
	if [ ! -f "$DIR_Mapping/$metagenome.bai" ]
	then
	((i=i%num_processes)); ((i++==0)) && wait
	/workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/Mapping_profiling_scripts/script-mapping_parallel.sh $metagenome &
	fi
done

wait

# END mapping time
ELAPSED="007-Mapping $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> ${mainDIR}/time-00-overall.tsv
