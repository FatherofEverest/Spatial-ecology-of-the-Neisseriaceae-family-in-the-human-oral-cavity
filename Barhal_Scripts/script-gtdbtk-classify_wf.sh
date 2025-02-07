#!/bin/bash

# variables
projectID=P_0003_Neisseriaceae
mainDIR=/workspace/jmarkwelchlab/$projectID
GTDB_genomes=$mainDIR/29_GTDBTK/G_0615/Genomes
OUT=$mainDIR/29_GTDBTK/G_0615/classify_out
TMP_DIR=$mainDIR/29_GTDBTK/G_0615/tmp
CUSTOM=$mainDIR/29_GTDBTK/G_0615/custom_taxonomy_file.txt
TH=5

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


# END time contigsDB
ELAPSED="GTDB-tk classify_wf $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"  
echo "$ELAPSED" >> $mainDIR/29_GTDBTK/G_0615/time.tsv
