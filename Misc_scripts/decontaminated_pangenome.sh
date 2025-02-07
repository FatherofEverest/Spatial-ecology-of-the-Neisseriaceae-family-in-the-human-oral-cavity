#!/bin/bash

export BLASTDB_LMDB_MAP_SIZE=1000000000

projectID=P_0003_Neisseriaceae
mainDIR=/workspace/jmarkwelchlab
DIR_Pangenome=$mainDIR/$projectID/09_PANGENOME
PanProject=${projectID}_G_0615_decontaminated-pangenome
PanDir=$DIR_Pangenome/$PanProject
GENOMES_decontaminated=$mainDIR/$projectID/DATA/${projectID}-decontaminated-contig_paths.txt
GENOMES_decontaminated_pangenome_DB=$PanDir/${PanProject}-GENOMES.db
PAN_DIR=$PanDir/${PanProject}-RESULTS
PAN_DB=$PanDir/${PanProject}-RESULTS/${PanProject}-PAN.db
ANI_DIR=$PanDir/${PanProject}-RESULTS/ANI_RESULTS
layersADD_decontaminated=DATA/$projectID-decontaminated-add_info.items.txt
PAN_DES_decontaminated=DATA/description_decontaminated_pangenome.txt
TH=10


# external Genome storage
#anvi-gen-genomes-storage -e $GENOMES_decontaminated -o $GENOMES_decontaminated_pangenome_DB

# external pangenome
anvi-pan-genome -g $GENOMES_decontaminated_pangenome_DB --I-know-this-is-not-a-good-idea --use-ncbi-blast --align-with muscle  --minbit 0.5 --mcl-inflation 10 -n ${PanProject} -o $PAN_DIR --num-threads $TH --enforce-hierarchical-clustering --description $PAN_DES_decontaminated

# add LAYERS
anvi-import-misc-data -t layers -p $PAN_DB $layersADD_decontaminated

# external genome similarity
anvi-compute-genome-similarity -e $GENOMES_decontaminated -o $ANI_DIR -p $PAN_DB --program pyANI --method ANIb --num-threads $TH --log-file initial_pangenome_pyANIlog


