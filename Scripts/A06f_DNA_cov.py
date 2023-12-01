
# A06f_DNA_cov.py ==============================================================

# setup ------------------—------------------—----------------------------------

import glob
import pandas as pd
import numpy as np

import os
filepath_wellmetadat = os.environ['metadat_well']
metadata_well = pd.read_csv(filepath_wellmetadat)

target_chroms = ["chr" + str(i) for i in range(1, 22)]
total_autosomal_bases = \
    pd.read_csv(os.environ['ref_chromsizes'],
                sep = "\t", header = None, index_col = 0).loc[target_chroms, 1].sum()



# gather metadata ------------------—------------------—------------------------

# unique coverage levels 
metadata_well = pd.read_csv(filepath_wellmetadat)

def parse_coverage_unique(filepath):
    try:
        percent_coverage = \
            pd.read_csv(filepath, delimiter = "\t", header = None, index_col=0)
        percent_coverage = percent_coverage.loc[
            np.intersect1d(target_chroms, percent_coverage.index), 1].sum() / total_autosomal_bases
    except:
        percent_coverage = pd.NA
    return(percent_coverage)

list_unique = [parse_coverage_unique(file) for file in metadata_well['A04f_txt_covtot']]
df_unique = pd.DataFrame(list_unique,
                        index = metadata_well['wellprefix'])

del(list_unique)
df_unique.columns = ["CoveragePerc1x"]
df_unique.to_csv("Metadata/A06f_DNA_cov_percent1x.tsv", sep='\t')
del(df_unique)

# total coverage levels 

def parse_coverage_total(filepath):
    try:
        total_cov_by_chr = pd.read_csv(filepath, delimiter = "\s+", header = None, index_col=1)
        coverage_XdivY = total_cov_by_chr.loc['chrX', 0] / total_cov_by_chr.loc['chrY', 0]
    except:
        coverage_XdivY = pd.NA
    return(coverage_XdivY)

list_total = [parse_coverage_total(file) for file in metadata_well['A04f_txt_covnsites']]
df_total = pd.DataFrame(list_total,
             index = metadata_well['wellprefix'])

del(list_total)
df_total.columns = ["CoverageXdivY"]
df_total.to_csv("Metadata/A06f_DNA_cov_chrXdivY.tsv", sep='\t')
del(df_total)
