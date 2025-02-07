#!/bin/bash

WD=/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter
DIR_SummaryPROF=$WD/08_PROFILE_SUMMARY
GENE_LEVEL_DETECTION=$WD/23_GENE_LEVEL_DETECTION
projectID=P_0622_Haemophilus_Aggregatibacter

# copy required files to a new directory named 23_GENE_LEVEL_DETECTION

echo "1. BEGIN - copying gene detection files"
for bin in `cat $WD/bin_list.txt`
  do
  for site in TD PP BM
    do
      cp $DIR_SummaryPROF/${projectID}_${site}-profile/bin_by_bin/$bin/${bin}-gene_detection.txt $GENE_LEVEL_DETECTION/${site}/${bin}-gene_detection.txt
  done
done
echo "1. FINISHED - copying gene detection files"

echo "2. BEGIN - making COMBO 30 files"
for bin in `cat $WD/bin_list.txt`
  do 
  $WD/SCRIPTS/indivMetaCombiner_DETECTION.py $bin 30
done
echo "2. FINISHED - making COMBO 30 files"

echo "3. BEGIN - moving COMBO 30 gene detection files"
mv $GENE_LEVEL_DETECTION/*-gene_detection.txt-COMBO-30 $GENE_LEVEL_DETECTION/COMBO_30_DATA/
echo "3. FINISHED - moving COMBO 30 gene detection files"

