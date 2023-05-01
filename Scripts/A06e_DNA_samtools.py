
# A06e_DNA_samtools.py =========================================================

# setup ------------------—------------------—----------------------------------

import glob
import pandas as pd

import os
filepath_wellmetadat = os.environ['metadat_well']
metadata_well = pd.read_csv(filepath_wellmetadat)

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


metadata_well = pd.read_csv(filepath_wellmetadat)

# paired-end
list_samstats_pe = [parse_samstats(file) for file in metadata_well["A04e_txt_samstats_PE"]]
df_samstats_pe = pd.DataFrame(list_samstats_pe,
             index = metadata_well["wellprefix"])

del(list_samstats_pe)
df_samstats_pe.to_csv("Metadata/A06e_DNA_samstats_PE.tsv", sep='\t')
del(df_samstats_pe)


# single-end
list_samstats_se = [parse_samstats(file) for file in metadata_well["A04e_txt_samstats_SE"]]
df_samstats_se = pd.DataFrame(list_samstats_se,
                     index = metadata_well["wellprefix"]
                             ).drop(["InsertSizeAvg", "InsertSizeSD"], axis = 1)

del(list_samstats_se)
df_samstats_se.to_csv("Metadata/A06e_DNA_samstats_SE.tsv", sep='\t')
del(df_samstats_se)


