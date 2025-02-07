#!/bin/bash

AA_seqs="/workspace/jmarkwelchlab/P_0003_Neisseriaceae/HIA_genes/HIA_genes_amino_acid_sequences_filtered_queries_and_known.fasta"
AA_msa_out="/workspace/jmarkwelchlab/P_0003_Neisseriaceae/HIA_genes/HIA_genes_amino_acid_sequences_filtered_queries_and_known_msa.fasta"

muscle -in $AA_seqs -out $AA_msa_out -maxiters 2
