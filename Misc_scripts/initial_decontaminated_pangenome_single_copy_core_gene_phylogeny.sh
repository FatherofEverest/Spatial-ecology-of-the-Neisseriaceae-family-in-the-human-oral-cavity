#!/bin/bash


projectID=P_0003_Neisseriaceae
mainDIR=/workspace/jmarkwelchlab
DIR_Pangenome=$mainDIR/$projectID/09_PANGENOME
PanProject=${projectID}-decontaminated-pangenome
PanDir=$DIR_Pangenome/$PanProject
GENOMES_decontaminated_pangenome_DB=$PanDir/${PanProject}-GENOMES.db
PAN_DIR=$PanDir/${PanProject}-RESULTS
PAN_DB=$PanDir/${PanProject}-RESULTS/${PanProject}-PAN.db
N_genomes=473



anvi-get-sequences-for-gene-clusters -g $GENOMES_decontaminated_pangenome_DB \
                                     -p $PAN_DB \
                                     -o $PAN_DIR/single_copy_core_genes-fasta \
                                     --max-num-genes-from-each-genome 1 \
                                     --min-num-genomes-gene-cluster-occurs $N_genomes \
                                     --min-geometric-homogeneity-index 1 \
                                     --concatenate-gene-clusters 
                                     
                                     

iqtree -s $PAN_DIR/single_copy_core_genes-fasta -nt AUTO -m WAG -bb 1000                                      
                                     
