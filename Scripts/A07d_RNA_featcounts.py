
# A07d_RNA_featcounts.py =======================================================

# setup ------------------------------------------------------------------------

import os
import pandas as pd

filepath_wellmetadat = os.environ['metadat_well']
metadata_well = pd.read_csv(filepath_wellmetadat)
batchnums = pd.unique(metadata_well['platenum'])


def parse_featurecounts(filepath):

    featc_summary = pd.read_csv(filepath, delimiter='\t')
    names_samples = [filename.split("/")[1] for filename in featc_summary.columns[1:]]
    names_features = featc_summary.iloc[ :, 0]

    # calc total read, tidy column names
    featc_summary = featc_summary.iloc[:, 1:].transpose()
    featc_summary = featc_summary.set_axis(names_samples, axis = 0).set_axis(names_features, axis = 1)
    featc_summary['TotalReadsFiltered'] = featc_summary.sum(axis = 1) # from A06c .Aligned.bam --> .Final.bam

    # other unassigned features should be zero (non-mapped filtered out)
    featc_summary = featc_summary[
        ['TotalReadsFiltered', 'Assigned', 'Unassigned_NoFeatures', 'Unassigned_Ambiguity']]
    
    return(featc_summary)



# gather metadata --------------------------------------------------------------

# gene-level -------------------------------------------------------------------

print("\n\ngene-level quants...")
filelist = pd.Series(["featurecounts_gene/PE_" + str(i) + ".summary" for i in batchnums])
boolean_fileexists = [os.path.exists(f) for f in filelist]
list_fcgene_PE = [ parse_featurecounts(f) for f in filelist[boolean_fileexists] ]
df_fcgene_PE = pd.concat(list_fcgene_PE)

filelist = pd.Series(["featurecounts_gene/SE1_" + str(i) + ".summary" for i in batchnums])
boolean_fileexists = [os.path.exists(f) for f in filelist]
list_fcgene_SE1 = [ parse_featurecounts(f) for f in filelist[boolean_fileexists] ]
df_fcgene_SE1 = pd.concat(list_fcgene_SE1)

filelist = pd.Series(["featurecounts_gene/SE2_" + str(i) + ".summary" for i in batchnums])
boolean_fileexists = [os.path.exists(f) for f in filelist]
list_fcgene_SE2 = [ parse_featurecounts(f) for f in filelist[boolean_fileexists] ]
df_fcgene_SE2 = pd.concat(list_fcgene_SE2)

fcgene_joined = \
    pd.concat([df_fcgene_PE.add_prefix("PE_"),
               df_fcgene_SE1.add_prefix("SE1_"),
               df_fcgene_SE2.add_prefix("SE2_")], axis = 1
               )
fcgene_joined.index.names = ["wellprefix"]

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
print(fcgene_joined.isna().sum().to_string())

print("number of duplicated wells:")
ndupe = fcgene_joined.index.duplicated().sum()
print(ndupe)

# final export
print("exporting Metadata/A07d_RNA_featcounts_gene.tsv of shape: {}".format(*fcgene_joined.shape))
fcgene_joined.to_csv("Metadata/A07d_RNA_featcounts_gene.tsv", sep = '\t')


# exon-level -------------------------------------------------------------------

print("\n\nexon-level quants...")
filelist = pd.Series(["featurecounts_exon/PE_" + str(i) + ".summary" for i in batchnums])
boolean_fileexists = [os.path.exists(f) for f in filelist]
list_fcexon_PE = [ parse_featurecounts(f) for f in filelist[boolean_fileexists] ]
df_fcexon_PE = pd.concat(list_fcexon_PE)

filelist = pd.Series(["featurecounts_exon/SE1_" + str(i) + ".summary" for i in batchnums])
boolean_fileexists = [os.path.exists(f) for f in filelist]
list_fcexon_SE1 = [ parse_featurecounts(f) for f in filelist[boolean_fileexists] ]
df_fcexon_SE1 = pd.concat(list_fcexon_SE1)

filelist = pd.Series(["featurecounts_exon/SE2_" + str(i) + ".summary" for i in batchnums])
boolean_fileexists = [os.path.exists(f) for f in filelist]
list_fcexon_SE2 = [ parse_featurecounts(f) for f in filelist[boolean_fileexists] ]
df_fcexon_SE2 = pd.concat(list_fcexon_SE2)

fcexon_joined = \
    pd.concat([df_fcexon_PE.add_prefix("PE_"),
               df_fcexon_SE1.add_prefix("SE1_"),
               df_fcexon_SE2.add_prefix("SE2_")], axis = 1
               )
fcexon_joined.index.names = ["wellprefix"]

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
print(fcexon_joined.isna().sum().to_string())

print("number of duplicated wells:")
ndupe = fcexon_joined.index.duplicated().sum()
print(ndupe)

# final export
print("exporting Metadata/A07d_RNA_featcounts_exon.tsv of shape: {}".format(*fcexon_joined.shape))
fcexon_joined.to_csv("Metadata/A07d_RNA_featcounts_exon.tsv", sep = '\t')
