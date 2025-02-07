#!/bin/bash
projectID=P_0003_Neisseriaceae
mainDIR=/workspace/jmarkwelchlab
DIR_Contigs=$mainDIR/$projectID/03_GENOMES_EDITED
DIR_ContigsDB=$mainDIR/$projectID/04_CONTIGS_DB
DIR_Assemblies=$mainDIR/$projectID/02_ASSEMBLIES
DIR_Data=$mainDIR/$projectID/DATA
genomesID=$DIR_Data/id_genomes.txt

for genome in `cat $genomesID`
do
EDITED_GENOMES=$DIR_Assemblies/${genome}-RAW.fa
CONTIGS=$DIR_Contigs/${genome}.fa
REPORT=$DIR_Contigs/${genome}.report.tsv
echo "reformatting ${genome}"
anvi-script-reformat-fasta -o $CONTIGS --simplify-names --prefix ${genome} --seq-type NT -r $REPORT $EDITED_GENOMES 
done


for genome in `cat $genomesID`
do
CONTIGS=$DIR_Contigs/${genome}.fa
contigsDB=$DIR_ContigsDB/${genome}-contigs.db
numThreads=5
echo "Generating anvio contigs db for ${genome}"
anvi-gen-contigs-database -f $CONTIGS -n ${genome} -o $contigsDB -T $numThreads 
done


