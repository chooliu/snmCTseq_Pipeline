
# A07c_RNA_featcounts.py =======================================================

# setup ------------------—------------------—----------------------------------

import pandas as pd

import os
filepath_wellmetadat = os.environ['metadat_well']
metadata_well = pd.read_csv(filepath_wellmetadat)

def parse_featurecounts(filepath):

    featc_summary = pd.read_csv(filepath, delimiter='\t')
    names_samples = [filename.split("/")[1] for filename in featc_summary.columns[1:]]
    names_features = featc_summary.iloc[ :, 0]

    # calc total read, tidy column names
    featc_summary = featc_summary.iloc[:, 1:].transpose()
    featc_summary = featc_summary.set_axis(names_samples, axis = 0).set_axis(names_features, axis = 1)
    featc_summary['TotalReadsFiltered'] = featc_summary.sum(axis = 1) # from A05c .Aligned.bam --> .Final.bam

    # other unassigned features should be zero (non-mapped filtered out)
    featc_summary = featc_summary[
        ['TotalReadsFiltered', 'Assigned', 'Unassigned_NoFeatures', 'Unassigned_Ambiguity']]
    
    return(featc_summary)



# gather metadata ------------------—------------------—------------------------

batchnums=pd.unique(metadata_well['platenum'])

# gene-level
list_fcgene_PE = [ parse_featurecounts("featurecounts_gene/PE_" + str(i) + ".summary")
                for i in batchnums ] 
df_fcgene_PE = pd.concat(list_fcgene_PE)

list_fcgene_SE1 = [ parse_featurecounts("featurecounts_gene/SE1_" + str(i) + ".summary")
                for i in batchnums ] 
df_fcgene_SE1 = pd.concat(list_fcgene_SE1)

list_fcgene_SE2 = [ parse_featurecounts("featurecounts_gene/SE2_" + str(i) + ".summary")
                for i in batchnums ] 
df_fcgene_SE2 = pd.concat(list_fcgene_SE2)

fcgene_joined = \
    pd.concat([df_fcgene_PE.add_prefix("PE_"),
               df_fcgene_SE1.add_prefix("SE1_"),
               df_fcgene_SE2.add_prefix("SE2_")], axis = 1
               )
fcgene_joined.index.names = ["wellprefix"]

fcgene_joined.to_csv("Metadata/A07c_RNA_featcounts_gene.tsv", sep='\t')


# intron-level
list_fcexon_PE = [ parse_featurecounts("featurecounts_exon/PE_" + str(i) + ".summary")
                for i in batchnums ] 
df_fcexon_PE = pd.concat(list_fcexon_PE)

list_fcexon_SE1 = [ parse_featurecounts("featurecounts_exon/SE1_" + str(i) + ".summary")
                for i in batchnums ] 
df_fcexon_SE1 = pd.concat(list_fcexon_SE1)

list_fcexon_SE2 = [ parse_featurecounts("featurecounts_exon/SE2_" + str(i) + ".summary")
                for i in batchnums ] 
df_fcexon_SE2 = pd.concat(list_fcexon_SE2)

fcexon_joined = \
    pd.concat([df_fcexon_PE.add_prefix("PE_"),
               df_fcexon_SE1.add_prefix("SE1_"),
               df_fcexon_SE2.add_prefix("SE2_")], axis = 1
               )
fcexon_joined.index.names = ["wellprefix"]

fcexon_joined.to_csv("Metadata/A07c_RNA_featcounts_exon.tsv", sep='\t')

