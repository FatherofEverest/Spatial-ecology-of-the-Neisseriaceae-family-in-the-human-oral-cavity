#!/bin/bash

genome=$1
site=$2

# Merged profile directory 
PROFILE_DB=$DIR_MergedPROF/${genome}/${genome}_${site}/PROFILE.db

# Contigs directory
CONTIGS_DB=$DIR_ContigsDB/${genome}-contigs.db

# Variability profile output
VAR_OUT=$DIR_Variability/${genome}/${site}/${genome}_SNVs.txt

# START timer per site
SECONDS=0

# Add deafault collection to profile db
anvi-script-add-default-collection -c $CONTIGS_DB -p $PROFILE_DB

# Generate variability reports
anvi-gen-variability-profile -c $CONTIGS_DB \
                             -p $PROFILE_DB \
                             -C DEFAULT \
                             -b EVERYTHING \
                             --compute-gene-coverage-stats \
                             --include-contig-names \
                             --include-split-names \
                             -o $VAR_OUT                             
# END timer per site
ELAPSED="${genome} ${site} $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> $mainDIR/27_VARIABILITY/Intra_species_diversity/time-variability_profile_per_site.tsv
