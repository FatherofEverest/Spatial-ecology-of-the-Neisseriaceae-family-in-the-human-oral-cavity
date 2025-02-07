#!/bin/bash


projectID=$1
mainDIR=/workspace/jmarkwelchlab/$projectID
CONTIGS_DB=04_CONTIGS_DB/${projectID}-contigs.db
TH=15

cd $mainDIR

# START time pfams
SECONDS=0
# run pfams
anvi-run-pfams -c $CONTIGS_DB --num-threads $TH
# END time Pfams
ELAPSED="008-Pfams $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> ${mainDIR}/time-00-overall.tsv

# START time NCBI COGs
SECONDS=0
# annotate NCBI COGs
anvi-run-ncbi-cogs -c $CONTIGS_DB --num-threads $TH --search-with blastp
# END time NCBI COGs
ELAPSED="009-NCBI_COGs $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> ${mainDIR}/time-00-overall.tsv

# START time KOfams
SECONDS=0
# annotate KOfamss
anvi-run-kegg-kofams -c $CONTIGS_DB --num-threads $TH --hmmer-program hmmscan
# END time NCBI COGs
ELAPSED="010-KOfams $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> ${mainDIR}/time-00-overall.tsv
