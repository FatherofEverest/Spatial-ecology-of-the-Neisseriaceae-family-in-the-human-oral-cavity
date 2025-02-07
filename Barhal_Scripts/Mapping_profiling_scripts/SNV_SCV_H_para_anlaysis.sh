#!/bin/bash

genome1=$1
genome2=$2

mainDIR=/workspace/jmarkwelchlab
projectID=P_0622_Haemophilus_Aggregatibacter
SUPP_profile=$mainDIR/$projectID/07_MERGED_PROFILE/${projectID}_PP/PROFILE.db
TD_profile=$mainDIR/$projectID/07_MERGED_PROFILE/${projectID}_TD/PROFILE.db
contigs=$mainDIR/$projectID/04_CONTIGS_DB/${projectID}-contigs.db
DIR_var=$mainDIR/$projectID/27_VARIABILITY
samples1=$DIR_var/SAMPLES_H_para_SUPP_best_breadth_depth.txt
samples2=$DIR_var/SAMPLES_H_para_TD_best_breadth_depth.txt


# SUPP


anvi-gen-variability-profile -p $SUPP_profile \
  -c $contigs \
  -C Genomes \
  -b $genome1 \
  --samples-of-interest $samples1 \
  --quince-mode \
  --compute-gene-coverage-stats \
  --include-contig-names \
  --include-split-names \
  --compute-gene-coverage-stats  \
  --include-site-pnps \
  -o $DIR_var/SNVs/SUPP/${genome1}_SNV_variability_report_SUPP.txt


anvi-gen-variability-profile -p $SUPP_profile \
  --engine CDN \
  -c $contigs \
  -C Genomes \
  -b $genome1 \
  --samples-of-interest $samples1 \
  --quince-mode \
  --compute-gene-coverage-stats \
  --include-contig-names \
  --include-split-names \
  --compute-gene-coverage-stats  \
  --include-site-pnps \
  -o $DIR_var/SCVs/SUPP/${genome1}_SCV_variability_report_SUPP.txt



anvi-gen-variability-profile -p $TD_profile \
  -c $contigs \
  -C Genomes \
  -b $genome2 \
  --samples-of-interest $samples2 \
  --quince-mode \
  --compute-gene-coverage-stats \
  --include-contig-names \
  --include-split-names \
  --compute-gene-coverage-stats  \
  --include-site-pnps \
  -o $DIR_var/SNVs/TD/${genome2}_SNV_variability_report_TD.txt



anvi-gen-variability-profile -p $TD_profile \
  --engine CDN \
  -c $contigs \
  -C Genomes \
  -b $genome2 \
  --samples-of-interest $samples2 \
  --quince-mode \
  --compute-gene-coverage-stats \
  --include-contig-names \
  --include-split-names \
  --compute-gene-coverage-stats  \
  --include-site-pnps \
  -o $DIR_var/SCVs/TD/${genome2}_SCV_variability_report_TD.txt
