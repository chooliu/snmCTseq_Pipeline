
# A05d_DNA_global_mCfracs.py ==========================================================

# setup ------------------------------------------------------------------------

import os
import glob
import pandas as pd

filepath_wellmetadat = os.environ['metadat_well']
metadata_well = pd.read_csv(filepath_wellmetadat)


# gather metadata --------------------------------------------------------------

filelist=pd.Series([ "Metadata/A04d_mCfrac_" + str(i) + ".tsv"
                for i in pd.unique(metadata_well['batchnum']) ])
boolean_fileexists = [os.path.exists(f) for f in filelist]
list_mCfracs = [ pd.read_csv(f, delimiter="\t") for f in filelist[boolean_fileexists] ] 
df_mCfracs = pd.concat(list_mCfracs)
df_mCfracs = df_mCfracs.rename(columns = {"Well" : "wellprefix"})
df_mCfracs = df_mCfracs.set_index("wellprefix")

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
print(df_mCfracs.isna().sum().to_string())

print("number of duplicated wells:")
ndupe = df_mCfracs.index.duplicated().sum()
print(ndupe)

# final export
print("exporting Metadata/A05d_DNA_global_mCfracs.tsv of shape: {}".format(*df_mCfracs.shape))
df_mCfracs.to_csv("Metadata/A05d_DNA_global_mCfracs.tsv", sep = '\t')
print("\n\n")
