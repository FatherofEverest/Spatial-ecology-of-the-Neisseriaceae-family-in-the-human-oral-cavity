#!/bin/bash

projectID=P_0003_Neisseriaceae
mainDIR=/workspace/jmarkwelchlab
DIR_NCBI=$mainDIR/$projectID/01_NCBI_GENOMES
downloadNCBI=$mainDIR/$projectID/DATA/03_P_0003_Neisseriaceae-download_url.txt

for gca_assembly in `cat $downloadNCBI`
do
gcaFile=$( echo $gca_assembly | awk -F'/' '{print $NF}' | sed -e 's/\.gz//' )
if [ ! -f $DIR_NCBI/${gcaFile} ]
then
echo "File not found - Downloading: ${gcaFile} "
wget -q -P $DIR_NCBI $gca_assembly
gunzip $DIR_NCBI/${gcaFile}.gz
echo "unzipped: ${gcaFile} "
fi
done
