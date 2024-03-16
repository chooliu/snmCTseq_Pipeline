
# A07f_RNA_picard.py ===========================================================

# setup ------------------------------------------------------------------------

import os
import glob
import pandas as pd

filepath_wellmetadat = os.environ['metadat_well']
metadata_well = pd.read_csv(filepath_wellmetadat)

# read picard log files
def parse_picard_rna(filepath):
    data_dedupe = pd.read_csv(filepath, delimiter = "\t",
                     comment = "#", nrows = 1).transpose()[0]
    return(data_dedupe)



# gather metadata --------------------------------------------------------------


# paired-end -------------------------------------------------------------------

print("\n\nPE logs...")
filelist = metadata_well['A06e_txt_picard_PE']
boolean_fileexists = [os.path.exists(f) for f in filelist]
list_picardrna_PE = [parse_picard_rna(f) for f in filelist[boolean_fileexists]]
df_picardrna_PE = pd.DataFrame(list_picardrna_PE,
                               index = metadata_well['wellprefix'][boolean_fileexists]
                               ).drop(["SAMPLE", "LIBRARY", "READ_GROUP"], axis = 1
                                      ).add_prefix("picard_")
df_picardrna_PE.columns = df_picardrna_PE.columns.str.lower()


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
print(df_picardrna_PE.isna().sum().to_string())

print("number of duplicated wells:")
ndupe = df_picardrna_PE.index.duplicated().sum()
print(ndupe)

# final export
print("exporting Metadata/A07f_RNA_picard_PE.tsv of shape: {}".format(*df_picardrna_PE.shape))
df_picardrna_PE.to_csv("Metadata/A07f_RNA_picard_PE.tsv", sep = '\t')
print("\n\n")





# single-end, read 1 -----------------------------------------------------------

print("\n\nSE1 logs...")
filelist = metadata_well['A06e_txt_picard_SE1']
boolean_fileexists = [os.path.exists(f) for f in filelist]
list_picardrna_SE1 = [parse_picard_rna(f) for f in filelist[boolean_fileexists]]
df_picardrna_SE1 = pd.DataFrame(list_picardrna_SE1,
                               index = metadata_well['wellprefix'][boolean_fileexists]
                               ).drop(["SAMPLE", "LIBRARY", "READ_GROUP"], axis = 1
                                      ).add_prefix("picard_")
df_picardrna_SE1.columns = df_picardrna_SE1.columns.str.lower()

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
print(df_picardrna_SE1.isna().sum().to_string())

print("number of duplicated wells:")
ndupe = df_picardrna_SE1.index.duplicated().sum()
print(ndupe)

# final export
print("exporting Metadata/A07f_RNA_picard_SE1.tsv of shape: {}".format(*df_picardrna_SE1.shape))
df_picardrna_SE1.to_csv("Metadata/A07f_RNA_picard_SE1.tsv", sep = '\t')
print("\n\n")



# single-end, read 2 -----------------------------------------------------------

print("\n\nSE2 logs...")
filelist = metadata_well['A06e_txt_picard_SE2']
boolean_fileexists = [os.path.exists(f) for f in filelist]
list_picardrna_SE2 = [parse_picard_rna(f) for f in filelist[boolean_fileexists]]
df_picardrna_SE2 = pd.DataFrame(list_picardrna_SE2,
                               index = metadata_well['wellprefix'][boolean_fileexists]
                               ).drop(["SAMPLE", "LIBRARY", "READ_GROUP"], axis = 1
                                      ).add_prefix("picard_")
df_picardrna_SE2.columns = df_picardrna_SE2.columns.str.lower()


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
print(df_picardrna_SE2.isna().sum().to_string())

print("number of duplicated wells:")
ndupe = df_picardrna_SE2.index.duplicated().sum()
print(ndupe)

# final export
print("exporting Metadata/A07f_RNA_picard_SE2.tsv of shape: {}".format(*df_picardrna_SE2.shape))
df_picardrna_SE2.to_csv("Metadata/A07f_RNA_picard_SE2.tsv", sep = '\t')
print("\n\n")
