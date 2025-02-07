#!/bin/bash

# variables
projectID=$1
genome_ID=$2
mainDIR=/workspace/jmarkwelchlab/$projectID
WD=$mainDIR/SNV_DIST_TEST
BAM_DIR=$mainDIR/27_VARIABILITY/Intra_species_diversity/05_MAPPING/${genome_ID}
CONTIGS_DB=$mainDIR/27_VARIABILITY/Intra_species_diversity/04_CONTIGS_DB/${genome_ID}-contigs.db
threadsRun=5
minContigSIZE=200


######### filtered bam files  #########
anvi-script-add-default-collection -c $CONTIGS_DB -p $WD/07_MERGED_PROFILE_FILTERED/PROFILE.db

# Generate variability report for filtered profiles
anvi-gen-variability-profile -c $CONTIGS_DB \
                             -p $WD/07_MERGED_PROFILE_FILTERED/PROFILE.db \
                             -C DEFAULT \
                             -b EVERYTHING \
                             --compute-gene-coverage-stats \
                             --include-contig-names \
                             --include-split-names \
                             -o $WD/08_VARIABILITY_FILTERED/variability.txt \
                             --quince-mode
                             
                             
                             
######### un-filtered bam files #########
# Add deafault collection to profile db
anvi-script-add-default-collection -c $CONTIGS_DB -p $WD/07_MERGED_PROFILE_UNFILTERED/PROFILE.db 

# Generate variability report for filtered profiles
anvi-gen-variability-profile -c $CONTIGS_DB \
                             -p $WD/07_MERGED_PROFILE_UNFILTERED/PROFILE.db \
                             -C DEFAULT \
                             -b EVERYTHING \
                             --compute-gene-coverage-stats \
                             --include-contig-names \
                             --include-split-names \
                             -o $WD/08_VARIABILITY_UNFILTERED/variability.txt \
                             --quince-mode
                                    
                   
                             
                             
                             
