
# A07d_RNA_samtools.py =========================================================

# setup ------------------—------------------—----------------------------------

import glob
import pandas as pd

import os
filepath_wellmetadat = os.environ['metadat_well']
metadata_well = pd.read_csv(filepath_wellmetadat)

# import samtools stats
def parse_samstats(filepath):

    term_dict = {
        'raw total sequences': f'FilteredSeqCount',
        'error rate': f'ErrorRate',
        'insert size average': f'InsertSizeAvg',
        'insert size standard deviation': f'InsertSizeSD',
        }

    with open(filepath) as report:
        report_dict = {}
        for line in report:
            try:
                lhs, rhs = line.split(':')
            except ValueError:
                continue
            try:
                report_dict[term_dict[lhs]] = rhs.strip().split('\t')[0]
            except KeyError:
                pass
            
    return(report_dict)




# gather metadata ------------------—------------------—------------------------

# paired-end
list_samstats_PE = [parse_samstats(file) for file in metadata_well["A05e_txt_samtools_PE"]]
df_samstats_PE = pd.DataFrame(list_samstats_PE,
             index = metadata_well["wellprefix"])
del(list_samstats_PE)
df_samstats_PE.to_csv("Metadata/A07d_RNA_samstats_PE.tsv", sep='\t')
del(df_samstats_PE)


# single-end, read 1
list_samstats_SE1 = [parse_samstats(file) for file in metadata_well["A05e_txt_samtools_SE1"]]
df_samstats_SE1 = pd.DataFrame(list_samstats_SE1,
             index = metadata_well["wellprefix"]
                               ).drop(["InsertSizeAvg", "InsertSizeSD"], axis = 1)
del(list_samstats_SE1)
df_samstats_SE1.to_csv("Metadata/A07d_RNA_samstats_SE1.tsv", sep='\t')
del(df_samstats_SE1)

# single-end, read 2
list_samstats_SE2 = [parse_samstats(file) for file in metadata_well["A05e_txt_samtools_SE2"]]
df_samstats_SE2 = pd.DataFrame(list_samstats_SE2,
             index = metadata_well["wellprefix"]
                              ).drop(["InsertSizeAvg", "InsertSizeSD"], axis = 1)
del(list_samstats_SE2)
df_samstats_SE2.to_csv("Metadata/A07d_RNA_samstats_SE2.tsv", sep='\t')
del(df_samstats_SE2)
