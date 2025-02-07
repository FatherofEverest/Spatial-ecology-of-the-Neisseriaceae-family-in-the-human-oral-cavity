#!/bin/bash

# variables
projectID=$1
genome=$2
G_ID=$3
mainDIR=/workspace/jmarkwelchlab/$projectID
PROJECT_CONTIGS=$mainDIR/03_GENOMES_EDITED/${G_ID}.fa
CONTIGS_DB=$mainDIR/27_VARIABILITY/Intra_species_diversity/04_CONTIGS_DB/${genome}-contigs.db
TH=4

cd $mainDIR

module load diamond

# START time contigsDB
SECONDS=0

# contigsDB
anvi-gen-contigs-database -f $PROJECT_CONTIGS -n $genome -o $CONTIGS_DB -T $TH

# END time contigsDB
ELAPSED="ContigsDB ${genome} $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"  
echo "$ELAPSED" >> $mainDIR/27_VARIABILITY/Intra_species_diversity/time-00-overall.tsv