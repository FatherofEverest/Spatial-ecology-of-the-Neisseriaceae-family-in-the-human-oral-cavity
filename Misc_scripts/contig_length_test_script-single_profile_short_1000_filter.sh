#!/bin/bash

sample=$1

# mappingVAR
BAM=$DIR_Mapping/$sample.bam

# single profle directory
SINGLE_PROFILE=$DIR_SinglePROF/$sample

# single profiling per sample
anvi-profile -i $BAM -c $CONTIGS_DB -o $SINGLE_PROFILE --min-contig-length 1000 -S $sample --num-threads $threadsRun --write-buffer-size-per-thread 30000

