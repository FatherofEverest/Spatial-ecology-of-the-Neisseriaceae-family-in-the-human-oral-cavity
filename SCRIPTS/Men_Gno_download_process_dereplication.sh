#!/bin/bash

# set up variables
projectID=P_0003_Neisseriaceae
mainDIR=/workspace/jmarkwelchlab
DIR_Data=$mainDIR/$projectID/DATA
DIR_NCBI=$mainDIR/$projectID/01_NCBI_GENOMES
DIR_NCBI_men_gno=$DIR_NCBI/MEN_GNO_GENOMES
DIR_Assemblies=$mainDIR/$projectID/02_ASSEMBLIES
DIR_Assemblies_men_gno=$mainDIR/$projectID/02_ASSEMBLIES/MEN_GNO_GENOMES
meningitidis_gonorrhoeae_metadata=$DIR_Data/filtered_meningitidis_gonorrhoeae_genbank_metadata.csv
meningitidis_gonorrhoeae_downloadNCBI=$DIR_Data/meningitidis_gonorrhoeae-download_url.txt
meningitidis_gonorrhoeae_projectMetadata=$DIR_Data/02_meningitidis_gonorrhoeae.csv
meningitidis_gonorrhoeae_rawNCBIid=$DIR_Data/meningitidis_gonorrhoeae-ncbi_raw_genomes_id.txt
meningitidis_gonorrhoeae_genomesNCBIid=$DIR_Data/meningitidis_gonorrhoeae-ncbi_genomes_id.txt
meningitidis_gonorrhoeae_genomesID=$DIR_Data/meningitidis_gonorrhoeae_id_genomes.txt
MEN_GNO_newHOMDid=$DIR_Data/meningitidis_gonorrhoeae-homd_id.txt
MEN_GNO_ITEMS=meningitidis_gonorrhoeae-add_info.items.txt
DIR_Contigs=$mainDIR/$projectID/03_GENOMES_EDITED
DIR_Contigs_men_gno=$DIR_Contigs/MEN_GNO_GENOMES
DIR_ContigsDB=$mainDIR/$projectID/04_CONTIGS_DB
DIR_ContigsDB_men_gno=$DIR_ContigsDB/MEN_GNO_GENOMES
MEN_GNO_GENOMES=$DIR_Data/MEN_GNO-contig_paths.txt
seqidHOMD=$mainDIR/HOMD_INFO/SEQID_info.txt
MEN_GNO_nameConversions=$DIR_Data/meningitidis_gonorrhoeae-name_conversions.txt
minContigSIZE=300


# make directories
mkdir $DIR_NCBI_men_gno $DIR_Assemblies_men_gno $DIR_Contigs_men_gno $DIR_ContigsDB_men_gno

# copy NCBI metadata for desired genomes (RefSeq)
cat $meningitidis_gonorrhoeae_metadata | awk 'BEGIN{FS=OFS=","}NR==1{print $0}NR>1{print $0}' > $meningitidis_gonorrhoeae_projectMetadata

# make list with gca path for download; GCA download link is in column 20
cat $meningitidis_gonorrhoeae_metadata | awk -F',' -v OFS="\t" 'NR>1{print $20}' | awk 'BEGIN{FS=OFS="/"}{print $0,$NF"_genomic.fna.gz"}' | sed -e 's/"//g' > $meningitidis_gonorrhoeae_downloadNCBI

#####  Download genomes ##### 
for gca_assembly in `cat $meningitidis_gonorrhoeae_downloadNCBI`
do
gcaFile=$( echo $gca_assembly | awk -F'/' '{print $NF}' | sed -e 's/\.gz//' )
if [ ! -f $DIR_NCBI_men_gno/${gcaFile} ]
then
echo "Downloading: ${gcaFile} "
wget -q -P $DIR_NCBI_men_gno $gca_assembly
gunzip $DIR_NCBI_men_gno/${gcaFile}.gz
echo "Successfully unzipped: ${gcaFile} "
fi
done


##### Rename genomes #####
# create a list of unique IDs containing a number and the assembly ID separated by "_id_"; use GCA assembly ID column 1
# (e.g. G_0001_id_GCA_000000000.0)
cat $meningitidis_gonorrhoeae_projectMetadata | awk -F',' -v OFS="\t" 'NR>1{ printf "G_""%04i_id_%s\n", NR-1,$1 }' | sed -e 's/"//g'  > $meningitidis_gonorrhoeae_rawNCBIid

# copy
cp $meningitidis_gonorrhoeae_rawNCBIid $meningitidis_gonorrhoeae_genomesNCBIid
cp $meningitidis_gonorrhoeae_rawNCBIid $meningitidis_gonorrhoeae_genomesID

# create new fasta files for each genome and 
# change names of original fasta file to number a number using the $genomesID file
# (e.g. GCA_901873365.1_Aggregatibacter_sp._BgEED05_genomic.fna -> G_0001-RAW.fa)
while IFS= read -r genomes_id
do
gca_id=$(echo $genomes_id | awk -F'_id_' '{print $2}')
old_name=$(find $DIR_NCBI_men_gno -name "$gca_id*")
NEW_NAME=$(echo "$genomes_id" | awk -F'_id_' -v new_dir="$DIR_Assemblies_men_gno" '{print new_dir"/"$1"-RAW.fa"}' )
cp $old_name $NEW_NAME
done < $meningitidis_gonorrhoeae_genomesID

# remove GCA number from $genomesID file to match new fasta file names
sed -i 's/_id_.*//' $meningitidis_gonorrhoeae_genomesID

#####  reformat genome deflines (edit this script to remove para) and build contigs db and make paths file
for genome in `cat $meningitidis_gonorrhoeae_genomesID`
do
EDITED_GENOMES=$DIR_Assemblies_men_gno/${genome}-RAW.fa
CONTIGS=$DIR_Contigs_men_gno/${genome}.fa
REPORT=$DIR_Contigs_men_gno/${genome}.report.tsv
anvi-script-reformat-fasta -l $minContigSIZE -o $CONTIGS --simplify-names --prefix ${genome} --seq-type NT -r $REPORT $EDITED_GENOMES 
done

# generate anvio contigs databases
for genome in `cat $meningitidis_gonorrhoeae_genomesID`
do
CONTIGS=$DIR_Contigs_men_gno/${genome}.fa
contigsDB=$DIR_ContigsDB_men_gno/${genome}-contigs.db
numThreads=10
anvi-gen-contigs-database -f $CONTIGS -n ${genome} -o $contigsDB -T $numThreads 
done

# get conversion names pangenome IDs to genome IDs
while IFS= read -r line
do
genomeNo=$( echo $line | awk -F'_id_' '{print $1}')
assemblyID=$( echo $line | awk -F'_id_' '{print $2}')
putGenoName=$( grep "$assemblyID" $meningitidis_gonorrhoeae_projectMetadata | awk -F',' -v OFS="\t" '{print $8,$9,$1}' | sed -e 's/"//g' | sed -e 's/ /_/g'  | sed -e 's/uncultured_//' | sed -e 's/\[//' | sed -e 's/\]//' | awk -F'\t' -v OFS="|" '{print $1,$2,$3}' | awk 'BEGIN{FS=OFS="_"}{ $1= substr($1,1,1)}1' | sed -e 's/_sp\./_sp/' | sed -E 's/(.*)\./\1POINT/' | sed -e "s/[^[:alnum:]|]/_/g" | sed -E 's/(.*)POINT/\1\./' | sed -e 's/_//' | sed -e 's/|/;/' |  sed -e 's/_.*;/;/' | sed -e 's/;/_str_/' | sed -e 's/|/_id_/' | sed -e 's/^.\{1\}/&_/' | sed -e 's/_str__id_/_str_NA_id_/' | sed -e 's/\./_/' | sed -e 's/strain_//g')
echo -e "$putGenoName\t$genomeNo" >> $MEN_GNO_nameConversions
done < $meningitidis_gonorrhoeae_rawNCBIid


# make meta data file for genome renaming
echo -e "item\tGenome_ID\tSpecies\tRefSeq\tAssembly_ID\tAssembly_ID_norm\tStrain\tGenome_in_HOMD\tHOMD_ID\tG_ID" > $MEN_GNO_ITEMS
for gcaID in `cat $MEN_GNO_nameConversions | cut -f1 | awk -F'_id_' '{print $2}' | rev | sed -e 's/_/./' | rev `
do
gcaID_ed=$(echo $gcaID | sed -e 's/\./_/')
itemsID=$(grep "$gcaID_ed" $MEN_GNO_nameConversions | cut -f1 | awk -F'_str_' -v OFS="\t" '{print $0,$0,$1}')
strainID=$( grep "$gcaID_ed" $MEN_GNO_nameConversions | cut -f1 | awk -F'_str_' -v OFS="\t" '{print $2}' | awk -F'_id_' -v OFS="\t" '{print $1}' )
refSeq=$(grep "$gcaID" $meningitidis_gonorrhoeae_projectMetadata | awk -F',' -v OFS="\t"  '{if ($22!="") print "Yes" ; else if ($16=="") print "No"}')
genomeHOMD=$( grep "$gcaID" $seqidHOMD)
if [ -z "$genomeHOMD" ]
then
inHOMD="no"
hmtID="NA"
else
inHOMD="yes"
hmtID=$( echo "$genomeHOMD"| awk 'BEGIN{FS=OFS="\t"}{print "HMT_"$2}')
fi
gID=$( grep "$gcaID_ed" $MEN_GNO_nameConversions | cut -f2 )
echo -e "$itemsID\t$refSeq\t$gcaID_ed\t$gcaID\t$strainID\t$inHOMD\t$hmtID\t$gID" >> $MEN_GNO_ITEMS
done


# Contig path file
echo -e 'name\tcontigs_db_path' > $MEN_GNO_GENOMES
for genome in `cat $meningitidis_gonorrhoeae_genomesID`
do
cat $MEN_GNO_nameConversions | grep "$genome" | awk -v workPath="$PWD" 'BEGIN{FS=OFS="\t"}{print $1,workPath"/04_CONTIGS_DB/MEN_GNO_GENOMES/"$2"-contigs.db"}' >> $MEN_GNO_GENOMES
done


# Re-do contig DBs that failed
anvi-gen-contigs-database -f $DIR_Contigs_men_gno/G_0175.fa -n G_0175 -o $DIR_ContigsDB_men_gno/G_0175-contigs.db -T 10 --force-overwrite         


for genome in `cat $meningitidis_gonorrhoeae_genomesID`
do
contigsDB=$DIR_ContigsDB_men_gno/${genome}-contigs.db
anvi-db-info $contigsDB | tee -a TMP_DB_INFO/Men_gen_contigs_db_info.txt
done


##### fastANI dereplication #####  

# 98% ANI
clusterize -n 10 -m jgiacomini@forsyth.org -l LOGS/Men_Gno_anvi-dereplicate-genomes_98_ANI.log anvi-dereplicate-genomes -e $MEN_GNO_GENOMES -o $DIR_Derep/meningitidis_gonorrhoeae_ANI_98_dereplication_fastANI --skip-fasta-report --program fastANI --similarity-threshold 0.98 --cluster-method simple_greedy --representative-method centrality --num-threads 20 --log-file $mainDIR/$projectID/LOGS/meningitidis_gonorrhoeae-anvi-dereplicate-fastANI_98.log

# 99% ANI
clusterize -n 10 -m jgiacomini@forsyth.org -l LOGS/Men_Gno_anvi-dereplicate-genomes_99_ANI.log anvi-dereplicate-genomes -e $MEN_GNO_GENOMES -o $DIR_Derep/meningitidis_gonorrhoeae_ANI_99_dereplication_fastANI --skip-fasta-report --program fastANI --similarity-threshold 0.99 --cluster-method simple_greedy --representative-method centrality --num-threads 10 --log-file $mainDIR/$projectID/LOGS/meningitidis_gonorrhoeae-anvi-dereplicate-fastANI_99.log






