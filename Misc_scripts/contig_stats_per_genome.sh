#!/bin/bash

projectID=P_0003_Neisseriaceae
mainDIR=/workspace/jmarkwelchlab
DIR_Contigs=$mainDIR/$projectID/04_CONTIGS_DB
DIR_OUT=$mainDIR/$projectID/CONTIG_LENGTH_TEST
DIR_Data=$mainDIR/$projectID/DATA
genomesID=$DIR_Data/id_genomes.txt


echo -e "Genome\tLongest_Contig\tShortest_Contig\tNum_Contigs_Greater_100kb\tNum_Contigs_Greater_50kb\tNum_Contigs_Greater_20kb\tNum_Contigs_Greater_10kb\tNum_Contigs_Greater_5kb\tNum_Contigs_Greater_2_5kb" > $DIR_OUT/Contig_stats.txt


for genome in `cat $genomesID`; do
  anvi-display-contigs-stats $DIR_Contigs/${genome}-contigs.db --report-as-text  -o $DIR_OUT/${genome}-contigsDB_stats.txt
  
  Longest_Contig=$(cat $DIR_OUT/${genome}-contigsDB_stats.txt | grep "Longest Contig" | cut -f2 )
  Shortest_Contig=$(cat $DIR_OUT/${genome}-contigsDB_stats.txt | grep "Shortest Contig" | cut -f2 )
  Num_Contigs_Greater_100kb=$(cat $DIR_OUT/${genome}-contigsDB_stats.txt | grep "Num Contigs > 100 kb" | cut -f2 )
  Num_Contigs_Greater_50kb=$(cat $DIR_OUT/${genome}-contigsDB_stats.txt | grep "Num Contigs > 50 kb" | cut -f2 )
  Num_Contigs_Greater_25kb=$(cat $DIR_OUT/${genome}-contigsDB_stats.txt | grep "Num Contigs > 25 kb" | cut -f2 )
  Num_Contigs_Greater_10kb=$(cat $DIR_OUT/${genome}-contigsDB_stats.txt | grep "Num Contigs > 10 kb" | cut -f2 )
  Num_Contigs_Greater_5kb=$(cat $DIR_OUT/${genome}-contigsDB_stats.txt | grep "Num Contigs > 5 kb" | cut -f2 )
  Num_Contigs_Greater_2_5kb=$(cat $DIR_OUT/${genome}-contigsDB_stats.txt | grep "Num Contigs > 2.5 kb" | cut -f2 )
  
  echo -e "$genome\t$Longest_Contig\t$Shortest_Contig\t$Num_Contigs_Greater_100kb\t$Num_Contigs_Greater_50kb\t$Num_Contigs_Greater_25kb\t$Num_Contigs_Greater_10kb\t$Num_Contigs_Greater_5kb\t$Num_Contigs_Greater_2_5kb" >> $DIR_OUT/Contig_stats.txt
  
  rm $DIR_OUT/${genome}-contigsDB_stats.txt

done


