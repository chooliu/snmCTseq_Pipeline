
# A07b_RNA_maprate.py ==========================================================

# setup ------------------—------------------—----------------------------------

import glob
import itertools
import re
import pandas as pd

import os
filepath_wellmetadat = os.environ['metadat_well']
metadata_well = pd.read_csv(filepath_wellmetadat)

def parse_star_report(filepath):

    """
    parse STAR.log output
    note that paired-end metrics usually fragments, versus reads
    """
    
    term_dict = {
        'Number of input reads': f'NumReadsIn',
        'Average input read length': f'AvgLengthIn',
        'Uniquely mapped reads number': f'NumReadsUniqueMapped',
        'Uniquely mapped reads %': f'PercentReadsUniqueMapped',
        'Average mapped length': f'AvgLengthMapped',
        'Number of splices: Total': f'NumTotSplices',
        'Number of splices: Annotated (sjdb)': f'NumAnnotSplices',
#         'Number of splices: GT/AG': f'NumGTAGSplices',
#         'Number of splices: GC/AG': f'NumGCAGSplices',
#         'Number of splices: AT/AC': f'NumATACSplices',
        'Mismatch rate per base, %': f'RateBaseMismatch',
        'Deletion rate per base': f'RateBaseDeletion',
        'Deletion average length': f'AvgLengthDeletion',
        'Insertion rate per base': f'RateBaseInsertion',
        'Insertion average length': f'AvgLengthInsertion',
#         'Number of reads mapped to multiple loci': f'NumReadsMultiMap',
        '% of reads mapped to multiple loci': f'PercentReadsMultiMap',
#         'Number of reads mapped to too many loci': f'NumReadsTooManyLoci',
        '% of reads mapped to too many loci': f'PercentReadsTooManyLoci',
#         'Number of reads unmapped: too many mismatches': f'NumReadsTooManyMismatch',
        '% of reads unmapped: too many mismatches':  f'PercentReadsTooManyMismatch',
#         'Number of reads unmapped: too short': f'NumReadsTooShort',
        '% of reads unmapped: too short': f'PercentReadsTooShort',
#         'Number of reads unmapped: other': f'NumReadsUnmappedOther',
        '% of reads unmapped: other': f'PercentReadsUnmappedOther',
#         'Number of chimeric reads': f'NumReadsChimeric',
#         '% of chimeric reads': f'PercentReadsChimeric',
    }
    
    with open(filepath) as report:
        report_dict = {}
        for line in report:
            try:
                lhs, rhs = line.split('|')
                lhs = lhs.strip()
            except ValueError:
                continue
            try:
                report_dict[term_dict[lhs]] = rhs.strip().strip('%')
            except KeyError:
                pass
            
    return(report_dict)
    
    


# gather metadata ------------------—------------------—------------------------

# paired-end
list_star_pe = [parse_star_report(file) for file in metadata_well['A05a_txt_star_PE']]
df_star_pe = pd.DataFrame(list_star_pe,
                          index= metadata_well['wellprefix'])
del(list_star_pe)
df_star_pe.to_csv("Metadata/A07b_RNA_maprate_PE.tsv", sep='\t')
del(df_star_pe)

# single-end, r1
list_star_SE1 = [parse_star_report(file) for file in metadata_well['A05a_txt_star_SE1']]
df_star_SE1 = pd.DataFrame(list_star_SE1,
                          index= metadata_well['wellprefix'])
del(list_star_SE1)
df_star_SE1.to_csv("Metadata/A07b_RNA_maprate_SE1.tsv", sep='\t')
del(df_star_SE1)

# single-end, r2
list_star_SE2 = [parse_star_report(file) for file in metadata_well['A05a_txt_star_SE2']]
df_star_SE2 = pd.DataFrame(list_star_SE2,
                          index= metadata_well['wellprefix'])
del(list_star_SE2)
df_star_SE2.to_csv("Metadata/A07b_RNA_maprate_SE2.tsv", sep='\t')
del(df_star_SE2)
