#!/bin/bash 


# input variables
projectID=$1
site=$2

# START of metapangenome time
SECONDS=0

# path to main directory
WORK_DIR=$mainDIR
# directory for metapangenome analysis
PAN_SITE=$DIR_Pangenome/HMP_${site}
# genomeDB for each site
GENOMES_STORAGE=$PAN_SITE/${projectID}-GENOMES.db
# pangenomeDB for each site
PAN_DB=$PAN_SITE/${projectID}-RESULTS/${projectID}-PAN.db
# summary directory for each site
panSummaryDIR=$DIR_SummaryPAN/HMP_${site}

# run summary for metapangenome
anvi-summarize -g $GENOMES_STORAGE -p $PAN_DB -C $panCOL -o $panSummaryDIR

# unzip summary file
gunzip $panSummaryDIR/${projectID}_gene_clusters_summary.txt.gz

# END of metapangenome time
ELAPSED="Summary_metapangenome_${site} $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> time-13-summary_metapangenome_per_site.tsv
