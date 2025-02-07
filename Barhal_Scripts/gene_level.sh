#!/bin/bash

projectID=P_0003_Neisseriaceae
mainDIR=/workspace/jmarkwelchlab
SUPP_profile=$mainDIR/$projectID/07_MERGED_PROFILE/${projectID}_PP/PROFILE.db
TD_profile=$mainDIR/$projectID/07_MERGED_PROFILE/${projectID}_TD/PROFILE.db
BM_profile=$mainDIR/$projectID/07_MERGED_PROFILE/${projectID}_KG/PROFILE.db
contigs=$mainDIR/$projectID/04_CONTIGS_DB/${projectID}-contigs.db
DIR_var=$mainDIR/$projectID/22_GENE_LEVEL
G_0213_pangenome_IDs=$mainDIR/$projectID/DATA/G_0213_pangenome_IDs.txt

for bin in `cat $G_0213_pangenome_IDs`
do 
anvi-script-gen-distribution-of-genes-in-a-bin -c $contigs -p $SUPP_profile -C Genomes -b $bin --fraction-of-median-coverage 0.25 
mv ${bin}-ENV-DETECTION.txt $DIR_var/SUPP/${bin}-ENV-DETECTION.txt
mv ${bin}-GENE-COVs.txt $DIR_var/SUPP/${bin}-GENE-COVs.txt
done

for bin in `cat $G_0213_pangenome_IDs`
do 
anvi-script-gen-distribution-of-genes-in-a-bin -c $contigs -p $TD_profile -C Genomes -b $bin --fraction-of-median-coverage 0.25 
mv ${bin}-ENV-DETECTION.txt $DIR_var/TD/${bin}-ENV-DETECTION.txt
mv ${bin}-GENE-COVs.txt $DIR_var/TD/${bin}-GENE-COVs.txt
done

for bin in `cat $G_0213_pangenome_IDs`
do 
anvi-script-gen-distribution-of-genes-in-a-bin -c $contigs -p $BM_profile -C Genomes -b $bin --fraction-of-median-coverage 0.25 
mv ${bin}-ENV-DETECTION.txt $DIR_var/KG/${bin}-ENV-DETECTION.txt
mv ${bin}-GENE-COVs.txt $DIR_var/KG/${bin}-GENE-COVs.txt
done

