#!/bin/bash
for file in  `cat /workspace/jmarkwelchlab/P_1221_Veillonella/samples_id-QC_IDs_836.txt`
do
samtools stats -X 05_MAPPING/${file}.bam 05_MAPPING/${file}.bai | grep ^SN | cut -f 2- > 05_MAPPING/${file}-STATS.txt
done
