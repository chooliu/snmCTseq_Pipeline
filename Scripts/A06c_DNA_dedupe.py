
# A06c_DNA_dedupe.py ===========================================================

# setup ------------------—------------------—----------------------------------

import glob
import pandas as pd
import numpy as np

import os
filepath_wellmetadat = os.environ['metadat_well']
metadata_well = pd.read_csv(filepath_wellmetadat)

# picard .log files
nulltable = np.array([pd.NA, pd.NA, pd.NA]) 

def parse_dedupe(filepath):
    try:
        data_dedupe = pd.read_csv(filepath, delimiter = "\t",
                         comment = "#", nrows = 1)[[
                             'UNPAIRED_READS_EXAMINED', 'READ_PAIRS_EXAMINED', 'PERCENT_DUPLICATION'
                         ]].transpose()[0]
        return(data_dedupe)
    except:
        return(nulltable)

tidy_name_dict = {'PERCENT_DUPLICATION' : 'picard_perc_dupe',
                  'READ_PAIRS_EXAMINED' : 'picard_npairsin',
                  'UNPAIRED_READS_EXAMINED' : 'picard_nreadsin'}



# gather metadata ------------------—------------------—------------------------

# paired end
list_picard_PE = [parse_dedupe(file) for file in metadata_well['A04a_txt_picard_PE']]
df_picard_PE = pd.DataFrame(list_picard_PE,
                            index = metadata_well['wellprefix']
                           ).rename(columns = tidy_name_dict
                           ).drop("picard_nreadsin", axis = 1)

del(list_picard_PE)
df_picard_PE.to_csv("Metadata/A06c_DNA_picard_PE.tsv", sep='\t')
del(df_picard_PE)

# single end
list_picard_SE = [parse_dedupe(file) for file in metadata_well['A04a_txt_picard_SE']]
df_picard_SE = pd.DataFrame(list_picard_SE,
                            index = metadata_well['wellprefix']
                           ).rename(columns = tidy_name_dict
                           ).drop("picard_npairsin", axis = 1)

del(list_picard_SE)
df_picard_SE.to_csv("Metadata/A06c_DNA_picard_SE.tsv", sep='\t')
del(df_picard_SE)



