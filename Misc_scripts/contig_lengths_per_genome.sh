#!/bin/bash

projectID=P_0003_Neisseriaceae
mainDIR=/workspace/jmarkwelchlab
DIR_Assemblies=$mainDIR/$projectID/02_ASSEMBLIES
DIR_OUT=$mainDIR/$projectID/CONTIG_LENGTH_TEST
DIR_Data=$mainDIR/$projectID/DATA
genomesID=$DIR_Data/id_genomes.txt


echo -e "genome\tcontig_length_threshold\tcontigs_removed\tTotal_no_contigs\tpercent_contigs_removed" > CONTIG_LENGTH_TEST/final_table.txt

for length in $(seq 100 100 1000); do 
   for genome in `cat $genomesID`; do
    anvi-script-reformat-fasta -l $length \
      -o $DIR_OUT/${genome}_${length}.fa \
      --simplify-names \
      --prefix $genome \
      --seq-type NT \
      $DIR_Assemblies/$genome-RAW.fa 2>&1 | tee $DIR_OUT/${genome}_${length}_terminal_output.txt 
    
    contigs_removed=$(cat $DIR_OUT/${genome}_${length}_terminal_output.txt | grep "Contigs removed" | sed -e 's/^.*: //g' | sed -e 's/ .*//g') 
    Total_no_contigs=$(cat $DIR_OUT/${genome}_${length}_terminal_output.txt | grep "Total num contigs" | sed -e 's/^.* //g') 
    percent_contigs_removed=$(cat $DIR_OUT/${genome}_${length}_terminal_output.txt | grep "Contigs removed" | sed -e 's/^.*(//g' | sed -e 's/ of all)//g') 

    echo -e "$genome\t$length\t$contigs_removed\t$Total_no_contigs\t$percent_contigs_removed" >> $DIR_OUT/final_table.txt 
    
    rm $DIR_OUT/${genome}_${length}.fa $DIR_OUT/${genome}_${length}_terminal_output.txt 
   done 
done


