#!/bin/bash

projectID=$1
site=$2

# Sample name
SAMPLE_NAME=${projectID}_${site}

# Merged profile directory
MERGED_PROFILE_DB=$DIR_MergedPROF/$SAMPLE_NAME/PROFILE.db

# SCG taxonomy matrix
SCG_TAX=DATA/${SAMPLE_NAME}-scg_taxonomy.tsv

# Summary directory
SUMMARY=$DIR_SummaryPROF/${SAMPLE_NAME}-profile

# START decompose profile time per site
SECONDS=0

anvi-import-collection -c $CONTIGS_DB -p $MERGED_PROFILE_DB -C $COLLECTION --contigs-mode $DECOMPOSE

# END decompose profile time  per sample
ELAPSED="${site} $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> ${mainDIR}/time-08-decompose_profile_per_site.tsv

# START summarize profile time per site
SECONDS=0

anvi-summarize -c $CONTIGS_DB -p $MERGED_PROFILE_DB --init-gene-coverages --report-aa-seqs-for-gene-calls -C $COLLECTION -o $SUMMARY

# END summarize profile time  per sample
ELAPSED="${site} $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> ${mainDIR}/time-10-summarized_profile_per_site.tsv
