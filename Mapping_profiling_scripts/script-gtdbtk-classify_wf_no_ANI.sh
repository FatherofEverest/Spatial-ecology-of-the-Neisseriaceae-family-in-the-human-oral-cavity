#!/bin/bash

# variables
projectID=P_0622_Haemophilus_Aggregatibacter
mainDIR=/workspace/jmarkwelchlab/$projectID
GTDB_genomes=$mainDIR/29_GTDBTK/Genomes
OUT=$mainDIR/29_GTDBTK/classify_out_no_ANI
CUSTOM=$mainDIR/29_GTDBTK/custom_taxonomy_file.txt
TH=1

cd $mainDIR

# START timer GTDB-tk
SECONDS=0

gtdbtk classify_wf \
    --genome_dir $GTDB_genomes \
    --out_dir $OUT \
    --skip_ani_screen \
    -x fa \
    --cpus $TH \
    --pplacer_cpus $TH


# END time contigsDB
ELAPSED="GTDB-tk classify_wf no ANI screen $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"  
echo "$ELAPSED" >> $mainDIR/29_GTDBTK/time.tsv
