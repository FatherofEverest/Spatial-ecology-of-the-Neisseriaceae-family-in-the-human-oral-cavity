#!/bin/bash

# Load variables
AA_seqs="/workspace/jmarkwelchlab/P_0003_Neisseriaceae/HIA_genes/HIA_genes_amino_acid_sequences_filtered_queries_121_and_known.fasta"
AA_msa_out="/workspace/jmarkwelchlab/P_0003_Neisseriaceae/HIA_genes/HIA_genes_amino_acid_sequences_filtered_queries_121_and_known_msa.fasta"

# Align sequences
muscle -in $AA_seqs -out $AA_msa_out

# Build ML tree with WAG model and bootstrap supoprt
iqtree -s $AA_msa_out -nt AUTO -m WAG -bb 1000 -o "WP_118780576.1" --redo-tree
