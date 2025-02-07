#!/bin/bash

WD=/workspace/jmarkwelchlab/P_0003_Neisseriaceae
DIR_SummaryPROF=$WD/08_PROFILE_SUMMARY
GENE_LEVEL_DETECTION=$WD/23_GENE_LEVEL_DETECTION
projectID=P_0003_Neisseriaceae

# copy required files to a new directory named 23_GENE_LEVEL_DETECTION
for bin in `cat $GENE_LEVEL_DETECTION/best_bins_IDs.txt`
  do
  for site in TD PP KG 
    do
      mkdir -p $GENE_LEVEL_DETECTION/${site}
      cp $DIR_SummaryPROF/${projectID}_${site}-profile/bin_by_bin/$bin/${bin}-gene_detection.txt $GENE_LEVEL_DETECTION/${site}/${bin}-gene_detection.txt
      echo "Finished copying data for ${bin} in ${site}"
  done
  
  # run bash script to make GENOME-gene_detection.txt-COMBO-30 files
  ./SCRIPTS/indivMetaCombiner_DETECTION.py ${bin} 30
  echo " Finished getting gene coverage data for ${bin}"

done
