
# A07c_RNA_dedupe.py ===========================================================

# setup ------------------------------------------------------------------------

import os
import glob
import pandas as pd
import numpy as np

filepath_wellmetadat = os.environ['metadat_well']
metadata_well = pd.read_csv(filepath_wellmetadat)

# picard .log files
nulltable = np.array([pd.NA, pd.NA, pd.NA]) 

def parse_picard_dedupe(filepath):
    try:
        data_dedupe = pd.read_csv(filepath, delimiter = "\t",
                         comment = "#", nrows = 1)[[
                             'UNPAIRED_READS_EXAMINED', 'READ_PAIRS_EXAMINED', 'PERCENT_DUPLICATION'
                         ]].transpose()[0]
        return(data_dedupe)
    except:
        print("error reading file: " + filepath)
        return(nulltable)

tidy_name_dict = {'PERCENT_DUPLICATION' : 'picard_perc_dupe',
                  'READ_PAIRS_EXAMINED' : 'picard_npairsin',
                  'UNPAIRED_READS_EXAMINED' : 'picard_nreadsin'}



# gather metadata ==============================================================

# paired-end -------------------------------------------------------------------

print("\n\nPE logs...")
filelist = metadata_well['A06b_bam_dedupe_PE']
boolean_fileexists = [os.path.exists(f) for f in filelist]
list_picard_PE = [parse_picard_dedupe(f) for f in filelist[boolean_fileexists]]
df_picard_PE = pd.DataFrame(list_picard_PE,
                        index = metadata_well['wellprefix'][boolean_fileexists])


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
print(df_picard_PE.isna().sum().to_string())

print("number of duplicated wells:")
ndupe = df_picard_PE.index.duplicated().sum()
print(ndupe)

# final export
print("exporting Metadata/A07c_RNA_picard_PE.tsv of shape: {}".format(*df_picard_PE.shape))
df_picard_PE.to_csv("Metadata/A07c_RNA_picard_PE.tsv", sep = '\t')
print("\n\n")


# single-end, read 1 -----------------------------------------------------------

print("\n\nSE1 logs...")
filelist = metadata_well['A06b_bam_dedupe_SE1']
boolean_fileexists = [os.path.exists(f) for f in filelist]
list_picard_SE1 = [parse_picard_dedupe(f) for f in filelist[boolean_fileexists]]
df_picard_SE1 = pd.DataFrame(list_picard_SE1,
                        index = metadata_well['wellprefix'][boolean_fileexists])

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
print(df_picard_SE1.isna().sum().to_string())

print("number of duplicated wells:")
ndupe = df_picard_SE1.index.duplicated().sum()
print(ndupe)

# final export
print("exporting Metadata/A07c_RNA_picard_SE1.tsv of shape: {}".format(*df_picard_SE1.shape))
df_picard_SE1.to_csv("Metadata/A07c_RNA_picard_SE1.tsv", sep = '\t')
print("\n\n")


# single-end, read 2 -----------------------------------------------------------

print("\n\nSE2 logs...")
filelist = metadata_well['A06b_bam_dedupe_SE2']
boolean_fileexists = [os.path.exists(f) for f in filelist]
list_picard_SE2 = [parse_picard_dedupe(f) for f in filelist[boolean_fileexists]]
df_picard_SE2 = pd.DataFrame(list_picard_SE2,
                        index = metadata_well['wellprefix'][boolean_fileexists])


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
print(df_picard_SE2.isna().sum().to_string())

print("number of duplicated wells:")
ndupe = df_picard_SE2.index.duplicated().sum()
print(ndupe)

# final export
print("exporting Metadata/A07c_RNA_picard_SE2.tsv of shape: {}".format(*df_picard_SE2.shape))
df_picard_SE2.to_csv("Metadata/A07c_RNA_picard_SE2.tsv", sep = '\t')
print("\n\n")
