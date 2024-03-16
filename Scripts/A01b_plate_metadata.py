


# ==============================================================================
# Scripts/A01b_plate_metadata.py
# should parse list of lane-merged plates -->
# extract plate-level metadata saved to $dir_proj/Metadata
# ==============================================================================

# recommend running interactively in python/Jupyter to check outputs,
# the relevant metadata parameters very likely to change between studies



# load packages ----------------------------------------------------------------

import glob
import sys
import pandas as pd
import os

# if running interactively, need to load some lines from snmCT_parameters.env
# or manually spec os.environ -- e.g., via os.environ['dir_proj'] = "mydirectory" or this below loop
# (check relative path of parameters.env file or change to absolute if below not working!)
envvar_needed = ['dir_proj', 'dir_originalfastq', 'metadat_plate']
try:
    os.environ['dir_proj']
except KeyError:
    envspec = pd.read_csv("../snmCT_parameters.env", sep = "=", comment="#", header = None
               ).set_axis(['varname', 'varpath'], axis = 1
               ).query('varname in @envvar_needed')
    for index, row in envspec.iterrows():
        os.environ[row["varname"]] = row["varpath"]
os.chdir(os.environ['dir_proj'])



# check fastq.gz names ---------------------------------------------------------

fastq_dir = os.environ['dir_originalfastq']
filepaths_raw_fastq = glob.glob(fastq_dir + "*fastq.gz")
print( filepaths_raw_fastq[0:4] )



# data.frame of plate names ----------------------------------------------------

# split before lane (L00[1-8]) to get unique plate names
plates_df = pd.DataFrame(
    {'plate' : pd.unique([filepath.split("/")[-1].split("_L")[0] for filepath in filepaths_raw_fastq])}
    ).sort_values('plate').reset_index(drop = True)

# .fastq name --> study specific metadata, may need customization # <--
# usually separated by "-"; example presented here is for IGVF cell lines
# if not custom may get "IndexError: list index out of range; ValueError: Transform function failed"
plates_df['datepool'] = plates_df['plate'].transform(lambda platename: platename.split("-")[0])
plates_df['sample'] = plates_df['plate'].transform(lambda platename: platename.split("-")[1])
plates_df['sort'] = plates_df['plate'].transform(lambda platename: platename.split("-")[2])
plates_df['plateindex'] = plates_df['plate'].transform(lambda platename: platename.split("-")[3])

plates_df['line'] = plates_df['sample'].transform(lambda platename: platename.split("D")[0])
plates_df['time'] = plates_df['sample'].transform(lambda platename: platename.split("D")[1])

# number each plate, "platenum" used for batch submission later on
# platenum indexed by 1-Nplates for compatibility with SGE (can't qsub -t 0)
plates_df['platenum'] = plates_df.index.astype(int) + 1
plates_df.index = plates_df.index.astype(int) + 1

# export to "Metadata/A01b_plate_metadata.csv" by default
print( plates_df.head() )
print ( plates_df.shape )
plates_df.to_csv(os.environ['metadat_plate'])
print("metadat_plate created.")


