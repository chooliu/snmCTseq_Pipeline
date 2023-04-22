
# A07e_RNA_picard.py ===========================================================

# setup ------------------—------------------—----------------------------------

import glob
import pandas as pd

import os
filepath_wellmetadat = os.environ['metadat_plate']
metadata_well = pd.read_csv(filepath_wellmetadat)

# read picard log files
def parse_picard_rna(filepath):
    data_dedupe = pd.read_csv(filepath, delimiter = "\t",
                     comment = "#", nrows = 1).transpose()[0]
    return(data_dedupe)



# gather metadata ------------------—------------------—------------------------

list_picard_PE = [parse_picard_rna(file) for file in metadata_well['A05e_txt_picard_PE']]
df_picard_PE = pd.DataFrame(list_picard_PE,
                            index = metadata_well['wellprefix']
                           ).drop(["SAMPLE", "LIBRARY", "READ_GROUP"], axis = 1
                           ).add_prefix("picard_")
df_picard_PE.columns = df_picard_PE.columns.str.lower()

del(list_picard_PE)
df_picard_PE.to_csv("Metadata/A07e_RNA_picard_PE.tsv", sep='\t')
del(df_picard_PE)

list_picard_SE1 = [parse_picard_rna(file) for file in metadata_well['A05e_txt_picard_SE1']]
df_picard_SE1 = pd.DataFrame(list_picard_SE1,
                            index = metadata_well['wellprefix']
                           ).drop(["SAMPLE", "LIBRARY", "READ_GROUP"], axis = 1
                           ).add_prefix("picard_")
df_picard_SE1.columns = df_picard_SE1.columns.str.lower()

del(list_picard_SE1)
df_picard_SE1.to_csv("Metadata/A07e_RNA_picard_SE1.tsv", sep='\t')
del(df_picard_SE1)

list_picard_SE2 = [parse_picard_rna(file) for file in metadata_well['A05e_txt_picard_SE2']]
df_picard_SE2 = pd.DataFrame(list_picard_SE2,
                            index = metadata_well['wellprefix']
                           ).drop(["SAMPLE", "LIBRARY", "READ_GROUP"], axis = 1
                           ).add_prefix("picard_")
df_picard_SE2.columns = df_picard_SE2.columns.str.lower()

del(list_picard_SE2)
df_picard_SE2.to_csv("Metadata/A07e_RNA_picard_SE2.tsv", sep='\t')
del(df_picard_SE2)
