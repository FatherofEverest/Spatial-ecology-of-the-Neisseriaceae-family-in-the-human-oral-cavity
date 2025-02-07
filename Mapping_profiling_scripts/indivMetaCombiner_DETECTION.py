#!/usr/bin/env python

from optparse import OptionParser
import pandas as pd
import os

parser = OptionParser()
(options, args) = parser.parse_args()

fileRoot=args[0]
keep=int(args[1])

# load detections
covs_td = pd.read_csv("23_GENE_LEVEL_DETECTION/TD/"+args[0]+"-gene_detection.txt", sep ="\t", index_col=False)
covs_bm = pd.read_csv("23_GENE_LEVEL_DETECTION/BM/"+args[0]+"-gene_detection.txt", sep ="\t", index_col=False)
covs_supp = pd.read_csv("23_GENE_LEVEL_DETECTION/PP/"+args[0]+"-gene_detection.txt", sep ="\t", index_col=False)



# sort by max median detection
covs_td = covs_td.reindex(covs_td.median().sort_values(ascending=False).index, axis=1)
covs_bm = covs_bm.reindex(covs_bm.median().sort_values(ascending=False).index, axis=1)
covs_supp = covs_supp.reindex(covs_supp.median().sort_values(ascending=False).index, axis=1)

# add in sites to samples
covs_td.columns.values[1:] = [x+"-TD" for x in covs_td.columns.values[1:]]
covs_bm.columns.values[1:] = [x+"-BM" for x in covs_bm.columns.values[1:]]
covs_supp.columns.values[1:] = [x+"-PP" for x in covs_supp.columns.values[1:]]

# make the new coverages
covs_keep = pd.concat([covs_td.iloc[:,0:(keep+1)], covs_bm.iloc[:,1:(keep+1)], covs_supp.iloc[:,1:(keep+1)]], axis=1)

# save it
covs_keep.to_csv("23_GENE_LEVEL_DETECTION/"+args[0]+"-gene_detection.txt-COMBO-"+str(keep), sep="\t", index=False)

