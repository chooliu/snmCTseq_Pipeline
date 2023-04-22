
# A06a_trimming.py =============================================================

# setup ------------------—------------------—----------------------------------

import re
import pandas as pd
import glob

import os
filepath_wellmetadat = os.environ['metadat_plate']
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



# gather metadata ------------------—------------------—------------------------

list_fastp = [parse_fastp_report(file) for file in metadata_well['A03a_json_fastp']]
df_fastp = pd.DataFrame(list_fastp,
                        index=metadata_well['wellprefix'])

del(list_fastp)
df_fastp.to_csv("Metadata/A06a_trimming.tsv", sep='\t')
del(df_fastp)
