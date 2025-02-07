#!/bin/bash

# variables
projectID=P_0622_Haemophilus_Aggregatibacter
mainDIR=/workspace/jmarkwelchlab/$projectID
GTDB_genomes=$mainDIR/29_GTDBTK/Genomes
OUT=$mainDIR/29_GTDBTK/classify_out
TMP_DIR=$mainDIR/29_GTDBTK/tmp
CUSTOM=$mainDIR/29_GTDBTK/custom_taxonomy_file.txt
TH=10

cd $mainDIR

# START timer GTDB-tk
SECONDS=0

gtdbtk classify_wf \
    --genome_dir $GTDB_genomes \
    --out_dir $OUT \
    --mash_db $OUT/gtdb-tk-r214.msh\
    -x fa \
    --cpus $TH \
    --pplacer_cpus $TH \
    --debug \
    --tmpdir $TMP_DIR

gtdbtk decorate \
    --input_tree $OUT/gtdbtk.bac120.classify.tree \
    --output_tree $OUT/Custom_gtdbtk.bac120.classify.tree \
    --custom_taxonomy_file $CUSTOM

# END time contigsDB
ELAPSED="GTDB-tk classify_wf $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"  
echo "$ELAPSED" >> $mainDIR/29_GTDBTK/time.tsv
