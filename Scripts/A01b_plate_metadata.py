
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

# if running interactively, check snmCT_parameters.env loaded or manually spec os.environ e.g.,
# os.environ['projdir'] ="/u/project/cluo/chliu/Analyses/IGVF"; os.chdir(os.environ['projdir'])
# os.environ['ref_dir'] = "/u/project/cluo/chliu/Genomes/human_gencode_v40"
# os.environ['dir_originalfastq'] = "/u/project/cluo/Shared_Datasets/IGVF/202208_Pilot/snmCT-seq/fastq/"
# os.environ['metadat_plate'] = "Metadata/A01b_plate_metadata.csv"



# check fastq.gz names ---------------------------------------------------------

fastq_dir = os.environ['dir_originalfastq']
filepaths_raw_fastq = glob.glob(fastq_dir + "*fastq.gz")
print( filepaths_raw_fastq[0:4] )


# data.frame of plate names ----------------------------------------------------

# split before lane (L00[1-4]) to get unique plate names
plates_df = pd.DataFrame(
    {'plate' : pd.unique([filepath.split("/")[-1].split("_L")[0] for filepath in filepaths_raw_fastq])}
    ).sort_values('plate').reindex()

# study specific metadata, usually separated by -
# example presented here is for IGVF cell lines
plates_df['dateseq'] = plates_df['plate'].transform(lambda platename: platename.split("-")[0])
plates_df['line'] = plates_df['plate'].transform(lambda platename: platename.split("-")[2])
plates_df['time'] = plates_df['plate'].transform(lambda platename: platename.split("-")[3])
plates_df['plateindex'] = plates_df['plate'].transform(lambda platename: platename.split("-")[4])

# number each plate, "platenum" used for batch submission later on
# platenum indexed by 1-Nplates for compatibility with SGE (can't qsub -t 0)
plates_df['platenum'] = plates_df.index.astype(int) + 1
plates_df.index = plates_df.index.astype(int) + 1

# export to "Metadata/A01b_plate_metadata.csv" by default
print( plates_df.head() )
print ( plates_df.shape )
plates_df.to_csv(os.environ['metadat_plate'])

