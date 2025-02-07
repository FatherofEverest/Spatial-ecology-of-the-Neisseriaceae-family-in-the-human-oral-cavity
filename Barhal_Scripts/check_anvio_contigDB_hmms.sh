#!/bin/bash
projectID=P_0003_Neisseriaceae
mainDIR=/workspace/jmarkwelchlab
DIR_ContigsDB=$mainDIR/$projectID/04_CONTIGS_DB
contigsDB=$DIR_ContigsDB

for contig in `ls $contigsDB`
  do
      anvi-display-contigs-stats  --report-as-text -o $mainDIR/$projectID/19_Contig_db_stats/$contig-stats  $DIR_ContigsDB/$contig
done

# check each contig-stats file 
echo -e "contig_db\tBacteria_71" > $mainDIR/$projectID/19_Contig_db_stats/Contig_bacteria_71_summary.txt

for contig in `ls $mainDIR/$projectID/19_Contig_db_stats`
  do
    echo "##### STARTING WITH $contig #####"
    cat 19_Contig_db_stats/$contig | grep 'contigs_db\|Bacteria_71' | sed 'n; n; d' > tmp
    contig_db=$(awk 'NR==1{print $2}' tmp)
    Bacteria_71=$(awk 'NR==2{print $2}' tmp)
    echo -e "$contig_db\t$Bacteria_71" >> $mainDIR/$projectID/19_Contig_db_stats/Contig_bacteria_71_summary.txt
    rm tmp
    echo "##### DONE WITH $contig #####"
done

