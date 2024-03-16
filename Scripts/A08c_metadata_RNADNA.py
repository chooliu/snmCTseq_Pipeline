
# A08c_metadata_RNADNA.py: compile A08a & A08b
# this final metadata file is what's typically used for QC

# setup ------------------------------------------------------------------------

import glob
import pandas as pd



# read and merge ---------------------------------------------------------------

def read_tbl_wrapper(filepath, prefix = ""):
    return(pd.read_csv(filepath, delimiter="\t", index_col = 0).add_prefix(prefix))

metadata_DNA = read_tbl_wrapper("Metadata/A08a_metadata_DNA.tsv", "DNA_")
metadata_RNA = read_tbl_wrapper("Metadata/A08b_metadata_RNA.tsv", "RNA_")

metadata_joined = \
    pd.concat([metadata_DNA,
               metadata_RNA.drop(
                   metadata_RNA.columns[metadata_RNA.columns.str.contains("Premap")],
                   axis = 1)],
              axis = 1)

# calculate few joint library size metrics
metadata_joined['Joint_TotalReadCount'] = \
    metadata_joined["DNA_Combined_Filtered_ReadCount"] + metadata_joined['RNA_Combined_Filtered_ReadCount']
metadata_joined['Joint_PercentRead_DNA'] = \
    metadata_joined["DNA_Combined_Filtered_ReadCount"] / metadata_joined['Joint_TotalReadCount']

metadata_joined['Joint_TotalFragmentCount'] = \
    metadata_joined["DNA_Combined_Filtered_FragmentCount"] + metadata_joined['RNA_Combined_Filtered_FragmentCount']
metadata_joined['Joint_PercentFragment_DNA'] = \
    metadata_joined["DNA_Combined_Filtered_FragmentCount"] / metadata_joined['Joint_TotalFragmentCount']

# RNA mapping rate, excluding DNA
metadata_joined['RNA_DNAAdj_ReadMappingRate'] = \
    metadata_joined['RNA_Combined_Filtered_ReadCount'] / \
    (metadata_joined['RNA_Combined_TotReadsIn'] - metadata_joined['DNA_Combined_Filtered_ReadCount'])

# DNA mapping rate, excluding RNA
metadata_joined['DNA_RNAAdj_ReadMappingRate'] = \
    metadata_joined['DNA_Combined_Filtered_ReadCount'] / \
    (metadata_joined['DNA_Combined_TotalReadsIn'] - metadata_joined['RNA_Combined_Filtered_ReadCount'])



# final DNA + RNA metadata -----------------------------------------------------

metadata_joined.to_csv("Metadata/A08c_metadata_RNADNA.tsv", sep = "\t")

