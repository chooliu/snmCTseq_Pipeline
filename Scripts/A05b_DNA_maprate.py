
# A05b_DNA_maprate.py ==========================================================



# setup ========================================================================

import os
import glob
import itertools
import re
import pandas as pd

filepath_wellmetadat = os.environ['metadat_well']
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





# gather metadata ==============================================================


# paired-end -------------------------------------------------------------------

print("\n\nPE logs...")
filelist = metadata_well['A04a_txt_bismark_PE']
boolean_fileexists = [os.path.exists(f) for f in filelist]
list_bismark_PE = [parse_bismark_report(f) for f in filelist[boolean_fileexists]]
df_bismark_PE = pd.DataFrame(list_bismark_PE,
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
print(df_bismark_PE.isna().sum().to_string())

print("number of duplicated wells:")
ndupe = df_bismark_PE.index.duplicated().sum()
print(ndupe)

# final export
print("exporting Metadata/A05b_DNA_maprate_PE.tsv of shape: {}".format(*df_bismark_PE.shape))
df_bismark_PE.to_csv("Metadata/A05b_DNA_maprate_PE.tsv", sep = '\t')
print("\n\n")





# read 1 singletons from trimming ----------------------------------------------

print("\n\nSE1trim logs...")
filelist = metadata_well['A04a_txt_bismark_SE1trim']
boolean_fileexists = [os.path.exists(f) for f in filelist]
list_bismark_SE1trim = [parse_bismark_report(f) for f in filelist[boolean_fileexists]]
df_bismark_SE1trim = pd.DataFrame(list_bismark_SE1trim,
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
print(df_bismark_SE1trim.isna().sum().to_string())

print("number of duplicated wells:")
ndupe = df_bismark_SE1trim.index.duplicated().sum()
print(ndupe)

# final export
print("exporting Metadata/A05b_DNA_maprate_SE1trim.tsv of shape: {}".format(*df_bismark_SE1trim.shape))
df_bismark_SE1trim.to_csv("Metadata/A05b_DNA_maprate_SE1trim.tsv", sep = '\t')
print("\n\n")





# read 2 singletons from trimming ----------------------------------------------

print("\n\nSE2trim logs...")
filelist = metadata_well['A04a_txt_bismark_SE2trim']
boolean_fileexists = [os.path.exists(f) for f in filelist]
list_bismark_SE2trim = [parse_bismark_report(f) for f in filelist[boolean_fileexists]]
df_bismark_SE2trim = pd.DataFrame(list_bismark_SE2trim,
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
print(df_bismark_SE2trim.isna().sum().to_string())

print("number of duplicated wells:")
ndupe = df_bismark_SE2trim.index.duplicated().sum()
print(ndupe)

# final export
print("exporting Metadata/A05b_DNA_maprate_SE2trim.tsv of shape: {}".format(*df_bismark_SE2trim.shape))
df_bismark_SE2trim.to_csv("Metadata/A05b_DNA_maprate_SE2trim.tsv", sep = '\t')
print("\n\n")





# read 1 singletons unampped in paired-end mode --------------------------------

print("\n\nSE1unmap logs...")
filelist = metadata_well['A04a_txt_bismark_SE1unmap']
boolean_fileexists = [os.path.exists(f) for f in filelist]
list_bismark_SE1unmap = [parse_bismark_report(f) for f in filelist[boolean_fileexists]]
df_bismark_SE1unmap = pd.DataFrame(list_bismark_SE1unmap,
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
print(df_bismark_SE1unmap.isna().sum().to_string())

print("number of duplicated wells:")
ndupe = df_bismark_SE1unmap.index.duplicated().sum()
print(ndupe)

# final export
print("exporting Metadata/A05b_DNA_maprate_SE1unmap.tsv of shape: {}".format(*df_bismark_SE1unmap.shape))
df_bismark_SE1unmap.to_csv("Metadata/A05b_DNA_maprate_SE1unmap.tsv", sep = '\t')
print("\n\n")



# read 2 singletons unampped in paired-end mode --------------------------------

print("\n\nSE2unmap logs...")
filelist = metadata_well['A04a_txt_bismark_SE2unmap']
boolean_fileexists = [os.path.exists(f) for f in filelist]
list_bismark_SE2unmap = [parse_bismark_report(f) for f in filelist[boolean_fileexists]]
df_bismark_SE2unmap = pd.DataFrame(list_bismark_SE2unmap,
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
print(df_bismark_SE2unmap.isna().sum().to_string())

print("number of duplicated wells:")
ndupe = df_bismark_SE2unmap.index.duplicated().sum()
print(ndupe)

# final export
print("exporting Metadata/A05b_DNA_maprate_SE2unmap.tsv of shape: {}".format(*df_bismark_SE2unmap.shape))
df_bismark_SE2unmap.to_csv("Metadata/A05b_DNA_maprate_SE2unmap.tsv", sep = '\t')
print("\n\n")
