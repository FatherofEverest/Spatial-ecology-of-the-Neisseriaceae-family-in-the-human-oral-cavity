#!/bin/bash

# variables
projectID=$1
site=$2

mainDIR=/workspace/jmarkwelchlab/$projectID
DIR_Contigs=$mainDIR/03_GENOMES_EDITED
DIR_ContigsDB=$mainDIR/04_CONTIGS_DB
DIR_MergedPROF=$mainDIR/07_MERGED_PROFILE
DIR_SummaryPROF=$mainDIR/08_PROFILE_SUMMARY
COLLECTION=Genomes
threadsRun=20
CONTIGS_DB=$DIR_ContigsDB/${projectID}-contigs.db

# Sample name
SAMPLE_NAME=${projectID}_${site}

# Merged profile directory
MERGED_PROFILE_DB=$DIR_MergedPROF/$SAMPLE_NAME/PROFILE.db

# Summary directory
SUMMARY=$DIR_SummaryPROF/${SAMPLE_NAME}-profile

# Summarize mapping results
anvi-summarize -c $CONTIGS_DB -p $MERGED_PROFILE_DB --init-gene-coverages --report-aa-seqs-for-gene-calls -C $COLLECTION -o $SUMMARY

