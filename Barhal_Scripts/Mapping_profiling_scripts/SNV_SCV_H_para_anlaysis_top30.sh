#!/bin/bash

genomeSUPP=$1
genomeTD=$2
genomeBM=$3
samplesSUPP=$4
samplesTD=$5
samplesBM=$6

mainDIR=/workspace/jmarkwelchlab
projectID=P_0622_Haemophilus_Aggregatibacter
SUPP_profile=$mainDIR/$projectID/07_MERGED_PROFILE/${projectID}_PP/PROFILE.db
TD_profile=$mainDIR/$projectID/07_MERGED_PROFILE/${projectID}_TD/PROFILE.db
BM_profile=$mainDIR/$projectID/07_MERGED_PROFILE/${projectID}_BM/PROFILE.db
contigs=$mainDIR/$projectID/04_CONTIGS_DB/${projectID}-contigs.db
DIR_var=$mainDIR/$projectID/27_VARIABILITY
samples=$DIR_var/DATA_melted_real_breadth_H_para_top30.txt


anvi-gen-variability-profile -p $SUPP_profile \
  -c $contigs \
  -C Genomes \
  -b $genomeSUPP \
  --samples-of-interest $samplesSUPP \
  --quince-mode \
  --compute-gene-coverage-stats \
  --include-contig-names \
  --include-split-names \
  --compute-gene-coverage-stats  \
  --include-site-pnps \
  -o $DIR_var/top30/SNVs/SUPP/${genome1}_SNV_variability_report_SUPP.txt


anvi-gen-variability-profile -p $TD_profile \
  -c $contigs \
  -C Genomes \
  -b $genomeTD \
  --samples-of-interest $samplesTD \
  --quince-mode \
  --compute-gene-coverage-stats \
  --include-contig-names \
  --include-split-names \
  --compute-gene-coverage-stats  \
  --include-site-pnps \
  -o $DIR_var/top30/SNVs/TD/${genome2}_SNV_variability_report_TD.txt


anvi-gen-variability-profile -p $BM_profile \
  -c $contigs \
  -C Genomes \
  -b $genomeBM \
  --samples-of-interest $samplesBM \
  --quince-mode \
  --compute-gene-coverage-stats \
  --include-contig-names \
  --include-split-names \
  --compute-gene-coverage-stats  \
  --include-site-pnps \
  -o $DIR_var/top30/SNVs/BM/${genome3}_SNV_variability_report_BM.txt
