#!/bin/bash

# build ML tree with WAG model and bootstrap supoprt
AA_msa_out="/workspace/jmarkwelchlab/P_0003_Neisseriaceae/HIA_genes/HIA_genes_amino_acid_sequences_filtered_queries_and_known_msa.fasta"

iqtree -s $AA_msa_out -nt AUTO -m WAG -bb 1000 -o "WP_118780576.1" --redo-tree
