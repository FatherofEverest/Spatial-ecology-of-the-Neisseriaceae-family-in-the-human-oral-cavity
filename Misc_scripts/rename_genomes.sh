#!/bin/bash

projectID=P_0003_Neisseriaceae
mainDIR=/workspace/jmarkwelchlab
DIR_NCBI=$mainDIR/$projectID/01_NCBI_GENOMES
DIR_Assemblies=$mainDIR/$projectID/02_ASSEMBLIES
DIR_Data=$mainDIR/$projectID/DATA
projectMetadata=$DIR_Data/02_${projectID}.csv
genomesNCBIid=$DIR_Data/08_P_0003_Neisseriaceae-ncbi_genomes_id.txt
rawNCBIid=$DIR_Data/04_P_0003_Neisseriaceae-ncbi_raw_genomes_id.txt
genomesID=$DIR_Data/id_genomes.txt

# create a list of unique IDs containing a number and the assembly ID separated by "_id_"; use GCA assembly ID column 9
# (e.g. G_0001_id_GCA_000000000.0)
cat $projectMetadata | awk -F',' -v OFS="\t" 'NR>1{ printf "G_""%04i_id_%s\n", NR-1,$12 }' | sed -e 's/"//g'  > $rawNCBIid

# copy
cp $rawNCBIid $genomesNCBIid
cp $rawNCBIid $genomesID

# create new fasta files for each genome and 
# change names of original fasta file to number a number using the $genomesID file
# (e.g. GCA_901873365.1_Aggregatibacter_sp._BgEED05_genomic.fna -> G_0001-RAW.fa)

while IFS= read -r genomes_id
do
gca_id=$(echo $genomes_id | awk -F'_id_' '{print $2}')
old_name=$(find $DIR_NCBI -name "$gca_id*")
NEW_NAME=$(echo "$genomes_id" | awk -F'_id_' -v new_dir="$DIR_Assemblies" '{print new_dir"/"$1"-RAW.fa"}' )
echo "copying ${old_name} to $NEW_NAME"
cp $old_name $NEW_NAME
done < $genomesID

# remove GCA number from $genomesID file to match new fasta file names
sed -i 's/_id_.*//' $genomesID


