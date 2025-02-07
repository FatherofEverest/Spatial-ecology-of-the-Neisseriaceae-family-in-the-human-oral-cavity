#!/bin/bash

# variables
projectID=P_0003_Neisseriaceae
mainDIR=/workspace/jmarkwelchlab
DIR_CLE=$mainDIR/$projectID/CONTIG_LENGTH_TEST
DIR_Contigs=$DIR_CLE/03_GENOMES_EDITED
DIR_Mapping=$DIR_CLE/05_MAPPING
EDITED_CONCAT_GENOMES_BOTH=$DIR_Contigs/Edited_Concatenated_both.fa
DIR_Reads=$mainDIR/P_0622_Haemophilus_Aggregatibacter/00_READS
sampleList=$DIR_CLE/Top_30_TD_samples.txt
threadsRun=10

# change to directory
cd $mainDIR

# generate REFERNCE-index for mapping
REFERENCE=$EDITED_CONCAT_GENOMES_BOTH
INDEX=$DIR_Mapping/Concatenated_reference_set

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
export DIR_CLE=$DIR_CLE

# START mapping time
SECONDS=0

num_processes=10

for metagenome in `cat $sampleList | awk 'NR>1{print $4}'`
do
        if [ ! -f "$DIR_Mapping/$metagenome.bai" ]
        then
        ((i=i%num_processes)); ((i++==0)) && wait
        /workspace/jmarkwelchlab/P_0003_Neisseriaceae/SCRIPTS/Contig_length_test_script-mapping_parallel.sh $metagenome &
        fi
done

wait

# END mapping time
ELAPSED="Mapping $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> ${DIR_CLE}/time.tsv

