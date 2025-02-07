#!/bin/bash

projectID=P_0003_Neisseriaceae
mainDIR=/workspace/jmarkwelchlab
DIR_ContigsDB=$mainDIR/$projectID/04_CONTIGS_DB
DIR_Data=$mainDIR/$projectID/DATA
genomesID=$DIR_Data/id_genomes.txt


module load diamond

for genome in `tail -n 8  $genomesID`
do
contigsDB=$DIR_ContigsDB/${genome}-contigs.db
numThreads=5
echo "Running anvio scg taxonomy for ${genome}"
anvi-run-scg-taxonomy -c $contigsDB -T $numThreads --write-buffer-size 30000
done

