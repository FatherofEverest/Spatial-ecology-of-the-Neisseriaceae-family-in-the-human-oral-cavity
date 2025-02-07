#!/bin/bash

projectID=P_0003_Neisseriaceae
mainDIR=/workspace/jmarkwelchlab
DIR_ContigsDB=$mainDIR/$projectID/04_CONTIGS_DB
DIR_Data=$mainDIR/$projectID/DATA
genomesID=$DIR_Data/id_genomes.txt

for i in `cat $genomesID` 
do
anvi-migrate $DIR_ContigsDB/${i}-contigs.db --migrate-safely
done

