#!/bin/bash


projectID=$1
mainDIR=/workspace/jmarkwelchlab/$projectID
PROJECT_CONTIGS=03_GENOMES_EDITED/${projectID}.fa
CONTIGS_DB=04_CONTIGS_DB/${projectID}-contigs.db
GENE_CALLS=11_GENE_CALLS/${projectID}-AAseq.fa
TH=15

cd $mainDIR
module load diamond

# START time scgTaxonomy
SECONDS=0
# get scgTaxonomy
anvi-run-scg-taxonomy -c $CONTIGS_DB --num-threads $TH --write-buffer-size 30000
# END time scgTaxonomyB
ELAPSED="005-scgTaxonomy $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> ${mainDIR}/time-00-overall.tsv

