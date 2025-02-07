#!/bin/bash

sample=$1

# START mapping time per sample
SECONDS=0

# readsVAR
READS=$DIR_Reads/$sample

# mappingVAR
MAP=$DIR_Mapping/$sample

# mapping Bowtie2 (very-sensitive, end-to-end, no-unal)
bowtie2 -x $INDEX -1 ${READS}_R1.fastq.gz -2 ${READS}_R2.fastq.gz -S $MAP.sam --end-to-end --very-sensitive --no-unal --threads $threadsRun

# SAM to BAM
samtools view -F 4 -bS $MAP.sam -o ${MAP}-RAW.bam --threads $threadsRun

# sort BAM
samtools sort -o $MAP.bam ${MAP}-RAW.bam --threads $threadsRun

# BAM index
samtools index $MAP.bam $MAP.bai

# remove intermediate files
rm ${MAP}-RAW.bam $MAP.sam

# END mapping time per sample
ELAPSED="${sample} $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> ${mainDIR}/time-05-mapping_per_sample.tsv
