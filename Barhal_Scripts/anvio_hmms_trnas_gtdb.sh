#!/bin/bash

projectID=P_0003_Neisseriaceae
mainDIR=/workspace/jmarkwelchlab
DIR_ContigsDB=$mainDIR/$projectID/04_CONTIGS_DB
DIR_Data=$mainDIR/$projectID/DATA
genomesID=$DIR_Data/id_genomes.txt

for genome in `cat $genomesID`
do
contigsDB=$DIR_ContigsDB/${genome}-contigs.db
numThreads=5
echo "Running anvio hmms for ${genome}"
anvi-run-hmms -c $contigsDB -T $numThreads 
done


# search for tRNAs in genomes (scan-trnas)
for genome in `cat $genomesID`
do
contigsDB=$DIR_ContigsDB/${genome}-contigs.db
numThreads=5
echo "Running anvio scan trnas for ${genome}"
anvi-scan-trnas -c $contigsDB -T $numThreads 
done

# GTDB taxonomy 

module load diamond
for genome in `cat $genomesID`
do
contigsDB=$DIR_ContigsDB/${genome}-contigs.db
numThreads=5
echo "Running anvio scg taxonomy for ${genome}"
anvi-run-scg-taxonomy -c $contigsDB -T $numThreads --write-buffer-size 30000 
done

