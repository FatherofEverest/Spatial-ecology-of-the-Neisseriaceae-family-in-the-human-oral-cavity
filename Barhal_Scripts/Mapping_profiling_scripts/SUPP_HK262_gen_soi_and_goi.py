#! /usr/bin/env python

import pandas as pd
import numpy as np
from anvio.terminal import Run

# -----------------------------------------------------------------
# Load in data
# -----------------------------------------------------------------

# Load up the gene coverage data
genes = pd.read_csv("08_PROFILE_SUMMARY/P_0622_Haemophilus_Aggregatibacter_PP-profile/bin_by_bin/H_parainfluenzae_str_HK262_id_GCA_000259485_1/H_parainfluenzae_str_HK262_id_GCA_000259485_1-gene_coverages.txt", sep='\t').set_index("gene_callers_id")

# Load up genome coverage data
genome = pd.read_csv("08_PROFILE_SUMMARY/P_0622_Haemophilus_Aggregatibacter_PP-profile/bins_across_samples/mean_coverage.txt", sep='\t')

# select only the bin of interest H_parainfluenzae_str_HK262_id_GCA_000259485_1
genome = genome[genome['bins'].str.contains("H_parainfluenzae_str_HK262_id_GCA_000259485_1")]  #.drop('bins', axis=1)

# Drop first column of dataframe
genome = genome.reset_index(drop=True)

# Subset metagenomes
metagenomes = genome[[col for col in genome.columns if col.startswith('PP_HC_HMP_')]]

# To determine samples and genes of interest, we do not include metagenomes

# remove column 1
del genome["bins"]

# transpose genome coverage data frame
genome = genome.T.rename(columns={0:'coverage'})

# -----------------------------------------------------------------
# Determine samples of interest
# -----------------------------------------------------------------

# Filter criteria: mean coverage should be > 10
samples_of_interest = genome.index[genome['coverage'] > 10].tolist()


# Additionally discard samples with coeff. of variation in gene coverages > 1.5
coeff_var = genes[samples_of_interest].std() / genes[samples_of_interest].mean()
samples_of_interest = coeff_var.index[coeff_var <= 1.25].tolist()

# Grab corresponding metagenomes for samples of interest
metagenomes_of_interest = ['PP_HC_HMP_' + sample for sample in samples_of_interest if 'PP_HC_HMP_' + sample in metagenomes]

# -----------------------------------------------------------------
# Determine genes of interest
# -----------------------------------------------------------------

# Establish lower and upper bound thresholds that all gene coverages must fall into
genome = genome[genome.index.isin(samples_of_interest)]
genome['greater_than'] = genome['coverage']/5
genome['less_than'] = genome['coverage']*50

genes = genes[samples_of_interest]

def is_core(row):
    return ((row > genome['greater_than']) & (row < genome['less_than'])).all()


genes_of_interest = genes.index[genes.apply(is_core, axis=1).astype(bool)].tolist()


# -----------------------------------------------------------------
# Summarize the info and write to file
# -----------------------------------------------------------------

run = Run()
run.warning("", header="Sample Info", nl_after=0, lc='green')
run.info('Num metagenomes', len(samples_of_interest))
run.info('metagenomes written to', 'SUPP_soi')

run.warning("", header="Gene Info", nl_after=0, lc='green')
run.info('Num genes', len(genes_of_interest))
run.info('genes written to', 'SUPP_goi', nl_after=1)

with open('SUPP_soi', 'w') as f:
    f.write('\n'.join([sample for sample in samples_of_interest]))
    f.write('\n')

with open('SUPP_goi', 'w') as f:
    f.write('\n'.join([str(gene) for gene in genes_of_interest]))
    f.write('\n')
