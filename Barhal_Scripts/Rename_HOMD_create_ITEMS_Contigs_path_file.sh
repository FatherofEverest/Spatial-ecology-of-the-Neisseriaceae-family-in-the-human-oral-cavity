#!/bin/bash

projectID=P_0003_Neisseriaceae
mainDIR=/workspace/jmarkwelchlab
DIR_Data=$mainDIR/$projectID/DATA
nameConversions=$DIR_Data/13_${projectID}-name_conversions.txt
rawNCBIid=$DIR_Data/new_04_${projectID}-ncbi_raw_genomes_id.txt
seqidHOMD=/workspace/jmarkwelchlab/HOMD_INFO/SEQID_info.txt
newHOMDid=$DIR_Data/07_${projectID}-homd_id.txt
ITEMS=${projectID}-add_info.items.txt
projectMetadata=$DIR_Data/02_${projectID}.csv
genomesID=$DIR_Data/id_genomes.txt

# get conversion names
# head is assembly_accession,host,isolation_source,Genus,Species,Strain_abrv,ftp_path
while IFS= read -r line
do
genomeNo=$( echo $line | awk -F'_id_' '{print $1}')
assemblyID=$( echo $line | awk -F'_id_' '{print $2}')
putGenoName=$( grep "$assemblyID" $projectMetadata | awk -F',' '{gsub(/ /, "_", $6); printf "%s_%s_str_%s_id_%s\n", substr($4,1,1), $5, $6, $1}' | sed -e 's/\./_/')
echo -e "$putGenoName\t$genomeNo" >> $nameConversions
done < $rawNCBIid

# find HOMD ID
for hypo_hmt in `cat $rawNCBIid | awk -F'_id_' '{print $2}' `
do
cat $seqidHOMD | grep "$hypo_hmt" | awk -v gca_id="$hypo_hmt" 'BEGIN{FS=OFS="\t"}{print $3"_"$4"_str_"$5"_id_"gca_id}' | sed -e's/ //g' | sed -e 's/_sp\./_sp_/' | awk 'BEGIN{FS=OFS="_"}{ $1= substr($1,1,1)}1' | awk '{if ($1 ~/HMT/) gsub("_sp","");print}' | awk -F'_id_' -v OFS="\t" '{print $NF,$0}' | sed -e 's/-/_/g' | grep 'HMT' | sed -e 's/HMT/sp_HMT_/' >> $newHOMDid
done

# edit name conversion with HOMD sps.
while IFS= read -r line
do
gcaID=$( echo $line | awk '{print $1}' | sed -e 's/\./_/' )
oldID=$( grep "$gcaID" $nameConversions | awk '{print $1}' )
newID=$( echo $line | awk '{print $2}' | sed -e 's/\./_/' )
sed -i "s/$oldID/$newID/" $nameConversions
done < $newHOMDid


# make meta data file for genome renaming
echo -e "item\tGenome_ID\tSpecies\tAssembly_ID\tAssembly_ID_norm\tStrain\tGenome_in_HOMD\tHOMD_ID\tG_ID" > $ITEMS
for gcaID in `cat $nameConversions | cut -f1 | awk -F'_id_' '{print $2}' | rev | sed -e 's/_/./' | rev `
do
gcaID_ed=$(echo $gcaID | sed -e 's/\./_/')
itemsID=$(grep "$gcaID_ed" $nameConversions | cut -f1 | awk -F'_str_' -v OFS="\t" '{print $0,$0,$1}')
strainID=$( grep "$gcaID_ed" $nameConversions | cut -f1 | awk -F'_str_' -v OFS="\t" '{print $2}' | awk -F'_id_' -v OFS="\t" '{print $1}' )
genomeHOMD=$( grep "$gcaID" $seqidHOMD)
if [ -z "$genomeHOMD" ]
then
inHOMD="no"
hmtID="NA"
else
inHOMD="yes"
hmtID=$( echo "$genomeHOMD"| awk 'BEGIN{FS=OFS="\t"}{print "HMT_"$2}')
fi
gID=$( grep "$gcaID_ed" $nameConversions | cut -f2 )
echo -e "$itemsID\t$gcaID_ed\t$gcaID\t$strainID\t$inHOMD\t$hmtID\t$gID" >> $ITEMS
done


# Contig path file
GENOMES=$mainDIR/$projectID/DATA/$projectID-contig_paths.txt
echo -e 'name\tcontigs_db_path' > $GENOMES
for genome in `cat $genomesID`
do
cat $nameConversions | grep "$genome" | awk -v workPath="$PWD" 'BEGIN{FS=OFS="\t"}{print $1,workPath"/04_CONTIGS_DB/"$2"-contigs.db"}' >> $GENOMES
done

