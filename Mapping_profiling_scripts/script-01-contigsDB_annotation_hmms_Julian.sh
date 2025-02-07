#!/bin/bash


projectID=$1
mainDIR=/workspace/jmarkwelchlab/$projectID
PROJECT_CONTIGS=03_GENOMES_EDITED/${projectID}.fa
CONTIGS_DB=04_CONTIGS_DB/${projectID}-contigs.db
GENE_CALLS=11_GENE_CALLS/${projectID}-AAseq.fa
TH=15

cd $mainDIR

module load diamond 

# START time contigsDB
SECONDS=0

# contigsDB
anvi-gen-contigs-database -f $PROJECT_CONTIGS -n $projectID -o $CONTIGS_DB -T $TH

# END time contigsDB
ELAPSED="001-contigsDB $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> ${mainDIR}/time-00-overall.tsv

# START time AAseq
SECONDS=0
# get AAseq
anvi-get-sequences-for-gene-calls -c $CONTIGS_DB --wrap 0 --get-aa-sequences -o $GENE_CALLS
# END time AAseq
ELAPSED="002-AAseq $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> ${mainDIR}/time-00-overall.tsv

# START time hmms
SECONDS=0
# run HMMs profiles
anvi-run-hmms -c $CONTIGS_DB --num-threads $TH
# END time hmms
ELAPSED="003-HMMs $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> ${mainDIR}/time-00-overall.tsv

# START time scan_tRNAs
SECONDS=0
# scan tRNAs
anvi-scan-trnas -c $CONTIGS_DB --num-threads $TH
# END time scan_tRNAs
ELAPSED="004-scan_tRNAs $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> ${mainDIR}/time-00-overall.tsv

# START time scgTaxonomy
SECONDS=0
# get scgTaxonomy
anvi-run-scg-taxonomy -c $CONTIGS_DB --num-threads $TH --write-buffer-size 30000
# END time scgTaxonomyB
ELAPSED="005-scgTaxonomy $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> ${mainDIR}/time-00-overall.tsv
