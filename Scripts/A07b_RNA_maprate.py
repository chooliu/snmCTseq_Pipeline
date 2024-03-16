
# A07b_RNA_maprate.py ==========================================================

# setup ------------------------------------------------------------------------

import os
import glob
import itertools
import re
import pandas as pd

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
    
    


# gather metadata --------------------------------------------------------------



# paired-end -------------------------------------------------------------------

print("\n\nPE logs...")
filelist = metadata_well['A06a_txt_star_PE']
boolean_fileexists = [os.path.exists(f) for f in filelist]
list_star_PE = [parse_star_report(f) for f in filelist[boolean_fileexists]]
df_star_PE = pd.DataFrame(list_star_PE,
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
print(df_star_PE.isna().sum().to_string())

print("number of duplicated wells:")
ndupe = df_star_PE.index.duplicated().sum()
print(ndupe)

# final export
print("exporting Metadata/A07b_RNA_maprate_PE.tsv of shape: {}".format(*df_star_PE.shape))
df_star_PE.to_csv("Metadata/A07b_RNA_maprate_PE.tsv", sep = '\t')
print("\n\n")



# single-end, r1 ---------------------------------------------------------------

print("\n\nSE1 logs...")
filelist = metadata_well['A06a_txt_star_SE1']
boolean_fileexists = [os.path.exists(f) for f in filelist]
list_star_SE1 = [parse_star_report(f) for f in filelist[boolean_fileexists]]
df_star_SE1 = pd.DataFrame(list_star_SE1,
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
print(df_star_SE1.isna().sum().to_string())

print("number of duplicated wells:")
ndupe = df_star_SE1.index.duplicated().sum()
print(ndupe)

# final export
print("exporting Metadata/A07b_RNA_maprate_SE1.tsv of shape: {}".format(*df_star_SE1.shape))
df_star_SE1.to_csv("Metadata/A07b_RNA_maprate_SE1.tsv", sep = '\t')
print("\n\n")



# single-end, r2 ---------------------------------------------------------------

print("\n\nSE2 logs...")
filelist = metadata_well['A06a_txt_star_SE2']
boolean_fileexists = [os.path.exists(f) for f in filelist]
list_star_SE2 = [parse_star_report(f) for f in filelist[boolean_fileexists]]
df_star_SE2 = pd.DataFrame(list_star_SE2,
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
print(df_star_SE2.isna().sum().to_string())

print("number of duplicated wells:")
ndupe = df_star_SE2.index.duplicated().sum()
print(ndupe)

# final export
print("exporting Metadata/A07b_RNA_maprate_SE2.tsv of shape: {}".format(*df_star_SE2.shape))
df_star_SE2.to_csv("Metadata/A07b_RNA_maprate_SE2.tsv", sep = '\t')
print("\n\n")
