#!/bin/bash

SUPP_profile=07_MERGED_PROFILE/P_0622_Haemophilus_Aggregatibacter_PP/PROFILE.db
TD_profile=07_MERGED_PROFILE/P_0622_Haemophilus_Aggregatibacter_TD/PROFILE.db
BM_profile=07_MERGED_PROFILE/P_0622_Haemophilus_Aggregatibacter_BM/PROFILE.db
contigs=04_CONTIGS_DB/P_0622_Haemophilus_Aggregatibacter-contigs.db
DIR_var=22_GENE_LEVEL
genomes_98ANI=DATA/id_genomes-98ANI.txt

#mkdir $DIR_var/SUPP $DIR_var/BM $DIR_var/TD

for bin in `cat $genomes_98ANI`
do 
anvi-script-gen-distribution-of-genes-in-a-bin -c $contigs -p $SUPP_profile -C Genomes -b $bin --fraction-of-median-coverage 0.25 
mv ${bin}-ENV-DETECTION.txt $DIR_var/SUPP/${bin}-ENV-DETECTION.txt
mv ${bin}-GENE-COVs.txt $DIR_var/SUPP/${bin}-GENE-COVs.txt
done


for bin in `cat $genomes_98ANI`
do 
anvi-script-gen-distribution-of-genes-in-a-bin -c $contigs -p $TD_profile -C Genomes -b $bin --fraction-of-median-coverage 0.25 
mv ${bin}-ENV-DETECTION.txt $DIR_var/TD/${bin}-ENV-DETECTION.txt
mv ${bin}-GENE-COVs.txt $DIR_var/TD/${bin}-GENE-COVs.txt
done

for bin in `cat $genomes_98ANI`
do 
anvi-script-gen-distribution-of-genes-in-a-bin -c $contigs -p $BM_profile -C Genomes -b $bin --fraction-of-median-coverage 0.25 
mv ${bin}-ENV-DETECTION.txt $DIR_var/BM/${bin}-ENV-DETECTION.txt
mv ${bin}-GENE-COVs.txt $DIR_var/BM/${bin}-GENE-COVs.txt
done

