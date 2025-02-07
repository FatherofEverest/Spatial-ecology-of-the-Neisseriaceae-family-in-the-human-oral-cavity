#!/bin/bash

export BLASTDB_LMDB_MAP_SIZE=1000000000

projectID=P_0003_Neisseriaceae
mainDIR=/workspace/jmarkwelchlab
DIR_Pangenome=$mainDIR/$projectID/09_PANGENOME
PanProject=${projectID}_G_0636-pangenome

PanDir=$DIR_Pangenome/$PanProject
GENOMES_G_0636=$mainDIR/$projectID/DATA/${projectID}-G_0636-contig_paths.txt
GENOMES_G_0636_pangenome_DB=$PanDir/${PanProject}-GENOMES.db
PAN_DIR=$PanDir/${PanProject}-RESULTS
PAN_DB=$PanDir/${PanProject}-RESULTS/${PanProject}-PAN.db
ANI_DIR=$PanDir/${PanProject}-RESULTS/ANI_RESULTS
PAN_DES_G_0636=$mainDIR/$projectID/DATA/description_G_0636_pangenome.txt
layersADD_G_0636=$mainDIR/$projectID/DATA/$projectID-G_0636-add_info.items.txt
TH=10


# external Genome storage
anvi-gen-genomes-storage -e $GENOMES_G_0636 -o $GENOMES_G_0636_pangenome_DB

# external pangenome
anvi-pan-genome -g $GENOMES_G_0636_pangenome_DB --I-know-this-is-not-a-good-idea --use-ncbi-blast --align-with muscle  --minbit 0.5 --mcl-inflation 10 -n ${PanProject} -o $PAN_DIR --num-threads $TH --enforce-hierarchical-clustering --description $PAN_DES_G_0636

# add LAYERS
anvi-import-misc-data -t layers -p $PAN_DB $layersADD_G_0636

# external genome similarity
anvi-compute-genome-similarity -e $GENOMES_G_0636 -o $ANI_DIR -p $PAN_DB --program fastANI --num-threads $TH --log-file G_0636_pangenome_fastANI_log


