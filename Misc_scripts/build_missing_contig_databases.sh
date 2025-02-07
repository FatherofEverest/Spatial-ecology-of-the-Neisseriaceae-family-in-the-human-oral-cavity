#!/bin/bash

projectID=P_0003_Neisseriaceae
mainDIR=/workspace/jmarkwelchlab
minContigSIZE=300
DIR_Contigs=$mainDIR/$projectID/03_GENOMES_EDITED
DIR_ContigsDB=$mainDIR/$projectID/04_CONTIGS_DB
DIR_Assemblies=$mainDIR/$projectID/02_ASSEMBLIES
DIR_Data=$mainDIR/$projectID/DATA
genomesID=$DIR_Data/id_genomes.txt

for genome in `tail -n8 $genomesID`
do
CONTIGS=$DIR_Contigs/${genome}.fa
contigsDB=$DIR_ContigsDB/${genome}-contigs.db
numThreads=5
anvi-gen-contigs-database -f $CONTIGS -n ${genome} -o $contigsDB -T $numThreads 
done

