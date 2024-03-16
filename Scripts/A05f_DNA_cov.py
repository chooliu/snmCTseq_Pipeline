
# A05f_DNA_cov.py ==============================================================

# setup ------------------------------------------------------------------------

import os
import glob
import pandas as pd
import numpy as np

filepath_wellmetadat = os.environ['metadat_well']
metadata_well = pd.read_csv(filepath_wellmetadat)

target_chroms = ["chr" + str(i) for i in range(1, 99)]
total_autosomal_bases = \
    pd.read_csv(os.environ['ref_chromsizes'],
                sep = "\t", header = None, index_col = 0)
total_autosomal_bases = \
    total_autosomal_bases.loc[np.intersect1d(target_chroms, total_autosomal_bases.index), 1].sum()


# gather metadata --------------------------------------------------------------


# extract autosomal ------------------------------------------------------------

target_chroms = ["chr" + str(i) for i in range(1, 99)]
autosomal_chroms = \
    pd.read_csv(os.environ['ref_chromsizes'],
                sep = "\t", header = None, index_col = 0)
autosomal_chroms = autosomal_chroms[autosomal_chroms.index.isin(target_chroms)]
total_autosomal_bases = autosomal_chroms.sum()
target_chroms = autosomal_chroms.index



# gather metadata: base-lvl unique coverage levels -----------------------------

print("processing autosomal num sites with at least 1-fold coverage.")
print("if any filenames printed below, potentially corrupt files:")
def parse_coverage_unique(filepath):
    try:
        percent_coverage = \
            pd.read_csv(filepath, delimiter = "\s+", header = None, index_col=1)
        percent_coverage = (
            percent_coverage.loc[
                np.intersect1d(target_chroms, percent_coverage.index), 0
            ].sum() / total_autosomal_bases).to_list()[0]
    except:
        print("'" + filepath + "'")
        percent_coverage = np.nan
    return(percent_coverage)

filelist = metadata_well['A04f_txt_covnsites']
boolean_fileexists = [os.path.exists(f) for f in filelist]
list_unique = [parse_coverage_unique(file) for file in filelist[boolean_fileexists]]
df_unique = pd.DataFrame(list_unique,
                        index = metadata_well['wellprefix'][boolean_fileexists])
df_unique.columns = ["CoveragePerc1x"]


# percent files missing
print("number of target files: " + str(len(filelist)))
print("fraction files missing: ")
print(round(1 - sum(boolean_fileexists)/len(boolean_fileexists), 3))
boolean_filemissing = [not f for f in boolean_fileexists]
if sum(boolean_filemissing) != 0:
    print("missing " + str(sum(boolean_filemissing)) + " files:")
    print(filelist[boolean_filemissing].to_string())

# column QC
print("number of NAs per column:")
print(df_unique.isna().sum().to_string())

print("number of duplicated wells:")
ndupe = df_unique.index.duplicated().sum()
print(ndupe)

# final export
print("exporting Metadata/A05f_DNA_cov_percent1x.tsv of shape: {}".format(*df_unique.shape))
df_unique.to_csv("Metadata/A05f_DNA_cov_percent1x.tsv", sep = '\t')
print("\n\n")



# total coverage levels for chrX/chrY ------------------------------------------

print("processing total coverage levels per chrom.")
print("if any filenames printed below, potentially corrupt files:")
def parse_coverage_total(filepath):
    try:
        total_cov_by_chr = pd.read_csv(filepath, delimiter = "\s+", header = None, index_col=0)
        if not(any(total_cov_by_chr.index=="chrX")) and (not any(total_cov_by_chr.index=="chrY")):
            coverage_XdivY = np.nan
        elif any(total_cov_by_chr.index=="chrX") and (not any(total_cov_by_chr.index=="chrY")):
            coverage_XdivY = np.inf
        else:
            coverage_XdivY = total_cov_by_chr.loc['chrX', ] / total_cov_by_chr.loc['chrY', ]
            coverage_XdivY = coverage_XdivY.tolist()[0]
    except:
        print("'" + filepath + "'")
        coverage_XdivY = np.nan
    return(coverage_XdivY)

filelist = metadata_well['A04f_txt_covtot']
boolean_fileexists = [os.path.exists(f) for f in filelist]
list_total = [parse_coverage_total(file) for file in filelist[boolean_fileexists]]
df_total = pd.DataFrame(list_total,
             index = metadata_well['wellprefix'][boolean_fileexists])
df_total.columns = ["CoverageXdivY"]


# percent files missing
print("number of target files: " + str(len(filelist)))
print("fraction files missing: ")
print(round(1 - sum(boolean_fileexists)/len(boolean_fileexists), 3))
boolean_filemissing = [not f for f in boolean_fileexists]
if sum(boolean_filemissing) != 0:
    print("missing " + str(sum(boolean_filemissing)) + " files:")
    print(filelist[boolean_filemissing].to_string())

# column QC
print("number of NAs per column:")
print(df_total.isna().sum().to_string())

print("number of duplicated wells:")
ndupe = df_total.index.duplicated().sum()
print(ndupe)

# final export
print("exporting Metadata/A05f_DNA_cov_chrXdivY.tsv of shape: {}".format(*df_total.shape))
df_total.to_csv("Metadata/A05f_DNA_cov_chrXdivY.tsv", sep = '\t')
print("\n\n")
