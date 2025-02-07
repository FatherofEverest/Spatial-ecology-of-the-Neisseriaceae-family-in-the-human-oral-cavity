#!/bin/bash

# input variables
projectID=$1
site=$2

# START of metapangenome time
SECONDS=0

# path to main directory
WORK_DIR=$mainDIR
# internal GENOMES.db
projectGenomesDB=$DIR_Pangenome/internal_annotated_${projectID}/${projectID}-GENOMES.db
# internal PAN.db
projectPanDIR=$DIR_Pangenome/internal_annotated_${projectID}/${projectID}-RESULTS
# profile collection
COLLECTION=Genomes
# contigsDB with path
CONTIGS_DB=$WORK_DIR/$DIR_ContigsDB/${projectID}-contigs.db
# profileDB with path
PROFILE_DB=$WORK_DIR/$DIR_MergedPROF/${projectID}_${site}/PROFILE.db
# directory for metapangenome analysis
PAN_SITE=$WORK_DIR/$DIR_Pangenome/HMP_${site}
# internal file for each site
INTERNAL=$PAN_SITE/internal_${projectID}_${site}.txt
# genomesDB for each site
GENOMES_STORAGE=$PAN_SITE/${projectID}-GENOMES.db
# pangenomeDB for each site
PAN_DB=$PAN_SITE/${projectID}-RESULTS/${projectID}-PAN.db

# make pangenome dir
mkdir $PAN_SITE

# copy internal GENOMES.db and PAN.db to site-pangenome directory
cp $projectGenomesDB $GENOMES_STORAGE
cp -r $projectPanDIR $PAN_SITE

# make internal file for each site
echo -e 'name\tbin_id\tcollection_id\tprofile_db_path\tcontigs_db_path' > $INTERNAL
cat $PAN_LAYERS | awk -v collection="$COLLECTION" -v profile_db="$PROFILE_DB" -v contigs_db="$CONTIGS_DB" 'BEGIN{FS=OFS="\t"}NR>1{print $1,$1,collection,profile_db,contigs_db}' >> $INTERNAL

# run metapangenome
anvi-meta-pan-genome -g $GENOMES_STORAGE --pan-db $PAN_DB -i $INTERNAL --fraction-of-median-coverage 0.25 --min-detection 0.5

# END of metapangenome time
ELAPSED="Metapangenome_${site} $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> time-12-metapangenome_per_site.tsv

