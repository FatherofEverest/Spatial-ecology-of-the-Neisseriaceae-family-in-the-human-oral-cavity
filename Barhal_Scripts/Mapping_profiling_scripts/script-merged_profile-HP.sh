#!/bin/bash

projectID=$1
site=$2

# Sample name
SAMPLE_NAME=${projectID}_${site}

# Merged profile directory
MERGED_PROFILE=$DIR_MergedPROF/$SAMPLE_NAME

# Single profiles per site
SINGLE_PROFILES_DB=$DIR_SinglePROF/${site}_*/PROFILE.db

# START merged profiling time per site
SECONDS=0

# merged profiling per sample
anvi-merge -c $CONTIGS_DB --enforce-hierarchical-clustering -o $MERGED_PROFILE -S $SAMPLE_NAME $SINGLE_PROFILES_DB

# END single profiling time  per sample
ELAPSED="${site} $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> ${mainDIR}/time-07-merged_profile_per_site.tsv
