#!/bin/bash


# make pangenome for one site (site to be edited later: SA = SV)
# this is the backbone pangenome, ani 
for site in SA
do
SECONDS=0
WORK_DIR=$PWD
TH=10
COLLECTION=Genomes
CONTIGS_DB=$WORK_DIR/$DIR_ContigsDB/${projectID}-contigs.db
PROFILE_DB=$WORK_DIR/$DIR_MergedPROF/PRACTICE/${projectID}_${site}/PROFILE.db
PAN_SITE=$WORK_DIR/$DIR_Pangenome/internal_annotated_${projectID}
INTERNAL=$PAN_SITE/internal_${projectID}.txt
GENOMES_STORAGE=$PAN_SITE/${projectID}-GENOMES.db
PAN_RESULTS=$PAN_SITE/${projectID}-RESULTS
PAN_DB=$PAN_RESULTS/${projectID}-PAN.db
ANIb_RESULTS=$PAN_SITE/ANIb-RESULTS
ANIt_RESULTS=$PAN_SITE/ANIt-RESULTS
mkdir $PAN_SITE
echo -e 'name\tbin_id\tcollection_id\tprofile_db_path\tcontigs_db_path' > $INTERNAL
cat $PAN_LAYERS | awk -v collection="$COLLECTION" -v profile_db="$PROFILE_DB" -v contigs_db="$CONTIGS_DB" 'BEGIN{FS=OFS="\t"}NR>1{print $1,$1,collection,profile_db,contigs_db}' >> $INTERNAL
anvi-gen-genomes-storage -i $INTERNAL --gene-caller prodigal -o $GENOMES_STORAGE
anvi-pan-genome -g $GENOMES_STORAGE --align-with muscle --use-ncbi-blast --minbit 0.5 --mcl-inflation 10 -n $projectID -o $PAN_RESULTS --num-threads $TH --enforce-hierarchical-clustering 
anvi-import-misc-data -p $PAN_DB -t layers $PAN_LAYERS
anvi-compute-genome-similarity -i $INTERNAL -o $ANIb_RESULTS --pan-db $PAN_DB --program pyANI --method ANIb --num-threads $TH
anvi-compute-genome-similarity -i $INTERNAL -o $ANIt_RESULTS --pan-db $PAN_DB --program pyANI --method TETRA --num-threads $TH



# hmms matrix
for hmms_source in Transfer_RNAs Ribosomal_RNA_28S Archaea_76 Ribosomal_RNA_16S Bacteria_71 Ptista_83 Ribosomal_RNA_18S Ribosomal_RNA_23S Ribosomal_RNA_5S Ribosomal_RNA_12S
do
HMMS_MATRIX=$DIR_Phylo/${projectID}-${hmms_source}-hits.tsv
anvi-script-gen-hmm-hits-matrix-across-genomes -i $INTERNAL --hmm-source ${hmms_source} -o $HMMS_MATRIX
done
ELAPSED="016-Pangenome_and_genome_similarity $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "$ELAPSED" >> time-11-internal_pangenome_per_site.tsv
echo "$ELAPSED" >> time-00-overall.tsv
done

