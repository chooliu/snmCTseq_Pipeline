
# A05e_DNA_samtools.py =========================================================

# setup ------------------------------------------------------------------------

import os
import glob
import pandas as pd

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




# gather metadata --------------------------------------------------------------

# paired-end -------------------------------------------------------------------

print("\n\npaired-end...")
filelist = metadata_well["A04e_txt_samstats_PE"]
boolean_fileexists = [os.path.exists(f) for f in filelist]
list_samstats_PE = [ parse_samstats(f) for f in filelist[boolean_fileexists] ] 
df_samstats_PE = pd.DataFrame(list_samstats_PE,
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
print(df_samstats_PE.isna().sum().to_string())

print("number of duplicated wells:")
ndupe = df_samstats_PE.index.duplicated().sum()
print(ndupe)

# final export
print("exporting Metadata/A05e_DNA_samstats_PE.tsv of shape: {}".format(*df_samstats_PE.shape))
df_samstats_PE.to_csv("Metadata/A05e_DNA_samstats_PE.tsv", sep = '\t')
print("\n\n")




# single-end -------------------------------------------------------------------

print("\n\nsingle-end...")
filelist = metadata_well["A04e_txt_samstats_SE"]
boolean_fileexists = [os.path.exists(f) for f in filelist]
list_samstats_SE = [ parse_samstats(f) for f in filelist[boolean_fileexists] ] 
df_samstats_SE = pd.DataFrame(list_samstats_SE,
                            index = metadata_well['wellprefix'][boolean_fileexists]
                               ).drop(['InsertSizeAvg', 'InsertSizeSD'], axis = 1)


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
print(df_samstats_SE.isna().sum().to_string())

print("number of duplicated wells:")
ndupe = df_samstats_SE.index.duplicated().sum()
print(ndupe)

# final export
print("exporting Metadata/A05e_DNA_samstats_SE.tsv of shape: {}".format(*df_samstats_SE.shape))
df_samstats_SE.to_csv("Metadata/A05e_DNA_samstats_SE.tsv", sep = '\t')
print("\n\n")
