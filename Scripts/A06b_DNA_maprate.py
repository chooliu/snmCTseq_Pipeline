
# A06b_DNA_maprate.py ==========================================================

# setup ------------------—------------------—----------------------------------

import glob
import itertools
import re
import pandas as pd

import os
filepath_wellmetadat = os.environ['metadat_plate']
metadata_well = pd.read_csv(filepath_wellmetadat)

    
def parse_bismark_report(filepath):

    """
    parse bismark.txt output
    adapted from YAP @ https://github.com/lhqing/cemba_data to include PE & SE output
    commented out term_dict lines of limited interest
    note that paired-end metrics usually yield fragments, versus reads
    """

    term_dict = {
        'Sequence pairs analysed in total': f'TotalReadPairsIn',
        'Sequences analysed in total': f'TotalReadsIn',
        'Number of paired-end alignments with a unique best hit': f'UniqueMappedPairs',
        'Number of alignments with a unique best hit from the different alignments': f'UniqueMappedReads',
        'Mapping efficiency': f'MappingRate',
#         'Sequence pairs with no alignments under any condition': f'UnmappedPairs',
#         'Sequences with no alignments under any condition': f'UnmappedReads',
#         'Sequences did not map uniquely': f'AmbigReads',
#         'Sequence pairs did not map uniquely': f'AmbigPairs',
#         'CT/GA/CT': f'ReadsOT',
#         'GA/CT/CT': f'ReadsOB',
#         'GA/CT/GA': f'ReadsCTOT',
#         'CT/GA/GA': f'ReadsCTOB',
#         'CT/CT': f'ReadsOT',
#         'CT/GA': f'ReadsOB',
#         'GA/CT': f'ReadsCTOT',
#         'GA/GA': f'ReadsCTOB',
#         'Total number of C\'s analysed': f'TotalC',
        'C methylated in CpG context': f'BismarkmCGRate',
        'C methylated in CHG context': f'BismarkmCHGRate',
        'C methylated in CHH context': f'BismarkmCHHRate',
        'C methylated in unknown context (CN or CHN)' : f'BismarkmCNCHNRate',
        'C methylated in Unknown context (CN or CHN)' : f'BismarkmCNCHNRate'
        }

    with open(filepath) as report:
        report_dict = {}
        for line in report:
            try:
                lhs, rhs = line.split(':')
            except ValueError:
                continue
            try:
                report_dict[term_dict[lhs]] = rhs.strip().split('\t')[0].strip('%')
            except KeyError:
                pass
            
    return(report_dict)



# gather metadata ------------------—------------------—------------------------

# paired-end bismark report file

list_bismark_PE = [parse_bismark_report(file) for file in metadata_well['A04a_txt_bismark_PE']]
df_bismark_PE = pd.DataFrame(list_bismark_PE,
                             index=metadata_well['wellprefix'])

del(list_bismark_PE)
df_bismark_PE.to_csv("Metadata/A06b_DNA_maprate_PE.tsv", sep='\t')
del(df_bismark_PE)




# read 1 singletons from trimming 

list_bismark_SE1trim = [parse_bismark_report(file) for file in metadata_well['A04a_txt_bismark_SE1trim']]
df_bismark_SE1trim = pd.DataFrame(list_bismark_SE1trim,
                             index=metadata_well['wellprefix'])

del(list_bismark_SE1trim)
df_bismark_SE1trim.to_csv("Metadata/A06b_DNA_maprate_SE1trim.tsv", sep='\t')
del(df_bismark_SE1trim)




# read 1 bismark unmapped in PE mode 

list_bismark_SE1unmap = [parse_bismark_report(file) for file in metadata_well['A04a_txt_bismark_SE1unmap']]
df_bismark_SE1unmap = pd.DataFrame(list_bismark_SE1unmap,
                             index=metadata_well['wellprefix'])

del(list_bismark_SE1unmap)
df_bismark_SE1unmap.to_csv("Metadata/A06b_DNA_maprate_SE1unmap.tsv", sep='\t')
del(df_bismark_SE1unmap)




# read 2 singletons from trimming 

list_bismark_SE2trim = [parse_bismark_report(file) for file in metadata_well['A04a_txt_bismark_SE2trim']]
df_bismark_SE2trim = pd.DataFrame(list_bismark_SE2trim,
                             index=metadata_well['wellprefix'])

del(list_bismark_SE2trim)
df_bismark_SE2trim.to_csv("Metadata/A06b_DNA_maprate_SE2trim.tsv", sep='\t')
del(df_bismark_SE2trim)



# read 2 bismark unmapped in PE mode 

list_bismark_SE2unmap = [parse_bismark_report(file) for file in metadata_well['A04a_txt_bismark_SE2unmap']]
df_bismark_SE2unmap = pd.DataFrame(list_bismark_SE2unmap,
                             index=metadata_well['wellprefix'])

del(list_bismark_SE2unmap)
df_bismark_SE2unmap.to_csv("Metadata/A06b_DNA_maprate_SE2unmap.tsv", sep='\t')
del(df_bismark_SE2unmap)



