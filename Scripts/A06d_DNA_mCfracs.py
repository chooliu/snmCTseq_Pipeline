
# A06d_DNA_mCfracs.py ==========================================================

# setup ------------------—------------------—----------------------------------

import glob
import pandas as pd

import os
filepath_wellmetadat = os.environ['metadat_plate']
metadata_well = pd.read_csv(filepath_wellmetadat)


# gather metadata ------------------—------------------—------------------------

list_mCfracs = [ pd.read_csv("Metadata/A04d_mCfrac_" + str(i) + ".tsv", delimiter="\t")
                for i in pd.unique(metadata_well['batchnum']) ] 
df_mCfracs = pd.concat(list_mCfracs)
df_mCfracs = df_mCfracs.rename(columns = {"Well" : "wellprefix"})

df_mCfracs.to_csv("Metadata/A06d_DNA_compiled_mCfracs.tsv", sep='\t', index = False)
