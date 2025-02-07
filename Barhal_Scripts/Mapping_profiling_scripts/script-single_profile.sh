#!/bin/bash

sample=$1

# START single profiling time per sample
SECONDS=0

# mappingVAR
BAM=$DIR_Mapping/$sample.bam

# single profle directory
SINGLE_PROFILE=$DIR_SinglePROF/$sample

# single profiling per sample
anvi-profile -i $BAM -c $CONTIGS_DB -o $SINGLE_PROFILE -S $sample --min-contig-length $minContigSIZE --profile-SCVs --num-threads $threadsRun --write-buffer-size-per-thread 30000

# END single profiling time  per sample
ELAPSED="${sample} $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> ${mainDIR}/time-06-single_profile_per_sample.tsv
