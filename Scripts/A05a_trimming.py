
# A05a_trimming.py =============================================================



# setup ========================================================================

import os
import re
import pandas as pd
import glob

filepath_wellmetadat = os.environ['metadat_well']
metadata_well = pd.read_csv(filepath_wellmetadat)

def parse_fastp_report(filepath):
    jsonfile = pd.read_json(filepath)
    dict_out = {
        'nreads_pretrim' : jsonfile['summary']['before_filtering']['total_reads'],
        'percreads_passtrim' : jsonfile['summary']['after_filtering']['total_reads'] /
              jsonfile['summary']['before_filtering']['total_reads'],
        'q20_pretrim' : jsonfile['summary']['before_filtering']['q30_rate'],
        'q20_posttrim' : jsonfile['summary']['after_filtering']['q30_rate'],
        'r1_len' : jsonfile['summary']['after_filtering']['read1_mean_length'],
        'r2_len' : jsonfile['summary']['after_filtering']['read2_mean_length'],
        'gc_perc' : jsonfile['summary']['after_filtering']['gc_content']}
    return(dict_out)



# gather metadata ==============================================================

print("\n\nfastp .json...")

filelist = metadata_well['A03a_json_fastp']
boolean_fileexists = [os.path.exists(f) for f in filelist]
list_fastp = [parse_fastp_report(f) for f in filelist[boolean_fileexists]]
df_fastp = pd.DataFrame(list_fastp,
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
print(df_fastp.isna().sum().to_string())

print("number of duplicated wells:")
ndupe = df_fastp.index.duplicated().sum()
print(ndupe)

# final export
print("exporting Metadata/A05a_trimming.tsv of shape: {}".format(*df_fastp.shape))
df_fastp.to_csv("Metadata/A05a_trimming.tsv", sep = '\t')
print("\n\n")
