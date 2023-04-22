
# A08b_final_metadat_RNA.py ====================================================
# assumes no changes to script output names from A06/A07

# setup ------------------—------------------—----------------------------------

import pandas as pd
import glob



# load tables ------------------—------------------—----------------------------

def read_tbl_wrapper(filepath, prefix = ""):
    return(pd.read_csv(filepath, delimiter="\t", index_col = 0).add_prefix(prefix))

df_trim_fastp = read_tbl_wrapper("Metadata/A06a_trimming.tsv", "Premap_")

df_star_pe = read_tbl_wrapper("Metadata/A07b_RNA_maprate_PE.tsv", "PE_")
df_star_se1 = read_tbl_wrapper("Metadata/A07b_RNA_maprate_SE1.tsv", "SE1_")
df_star_se2 = read_tbl_wrapper("Metadata/A07b_RNA_maprate_SE2.tsv", "SE2_")

df_featurecounts_gene = read_tbl_wrapper("Metadata/A07c_RNA_featcounts_gene.tsv", "Gene_")
df_featurecounts_exon = read_tbl_wrapper("Metadata/A07c_RNA_featcounts_exon.tsv", "Exon_")

df_samstat_pe = read_tbl_wrapper("Metadata/A07d_RNA_samstats_PE.tsv", "PE_")
df_samstat_se1 = read_tbl_wrapper("Metadata/A07d_RNA_samstats_SE1.tsv", "SE1_")
df_samstat_se2 = read_tbl_wrapper("Metadata/A07d_RNA_samstats_SE2.tsv", "SE2_")

df_picard_pe = read_tbl_wrapper("Metadata/A07e_RNA_picard_PE.tsv", "PE_")
df_picard_se1 = read_tbl_wrapper("Metadata/A07e_RNA_picard_SE1.tsv", "SE1_")
df_picard_se2 = read_tbl_wrapper("Metadata/A07e_RNA_picard_SE2.tsv", "SE2_")


metadata_rna = \
    pd.concat([df_trim_fastp,
               df_star_pe, df_star_se1, df_star_se2,
              df_featurecounts_gene, df_featurecounts_exon,
              df_samstat_pe, df_samstat_se1, df_samstat_se2,
              df_picard_pe, df_picard_se1, df_picard_se2],
              axis = 1)
metadata_rna = metadata_rna.apply(pd.to_numeric, errors='coerce')



# few 'combined' PE & SE metadata stats -------—--------------------------------
# some attempts at combined/weighted metrics

# notes: (i) most output above reports paired-end read-pairs/fragments versus reads, hence *2
#       (ii) reads that are attempted to be remapped in "stage-two" SE mapping
#            are subsets of "PE_TotalReadsIn"; thus can't simply sum(total mapped) / sum(reads in)
#     & (iii) "PE-mapping singletons" can also be generated where R1 assigned DNA while R2 assigned RNA,
#           or with STAR, R1 can map while R2 doesn't (e.g., -f 0x0048) and these are retained in PE.bam;
#           thus # based solely on the mapper's output miscount mapping singletons as PE instead of SE
#           & the single-end mapping & filtering rate interpretations are slightly off versus bismark
#           [may write QC scripts based on the _annotations files to make this more accurate in the future]


# estimated total mapping rate (at read level) ---------------------------------

# total reads in approximate because uses the rounded value of % reads mapped
# (& other caveats above); in practice better to use "Combined_TotalReadsIn" from mC calc in A07b

metadata_rna['Combined_TotReadsIn'] = \
    metadata_rna['PE_NumReadsIn'] * (metadata_rna['PE_PercentReadsUniqueMapped']/100) * 2 + \
    metadata_rna['SE1_NumReadsIn'] + \
    metadata_rna['SE2_NumReadsIn']

metadata_rna['Combined_NumReadsUniqueMapped'] = \
    metadata_rna['PE_NumReadsUniqueMapped']*2 + \
    metadata_rna['SE1_NumReadsUniqueMapped'] + \
    metadata_rna['SE2_NumReadsUniqueMapped']

metadata_rna['Combined_MappingRate'] = \
    metadata_rna['Combined_NumReadsUniqueMapped'] / metadata_rna['Combined_TotReadsIn']

# percent of mapped reads that pass filtering ----------------------------------
# note: SE1_TotalReadsFiltered should = SE1_FilteredSeqCount

metadata_rna['PE_Filtered_Rate'] = \
    metadata_rna['PE_FilteredSeqCount'] / \
    metadata_rna['PE_NumReadsUniqueMapped']
metadata_rna['SE1_Filtered_Rate'] = \
    metadata_rna['SE1_FilteredSeqCount'] / \
    metadata_rna['SE1_NumReadsUniqueMapped']
metadata_rna['SE2_Filtered_Rate'] = \
    metadata_rna['SE2_FilteredSeqCount'] / \
    metadata_rna['SE2_NumReadsUniqueMapped']

metadata_rna['Combined_Filtered_ReadCount'] = \
    metadata_rna['PE_FilteredSeqCount']*2 + \
    metadata_rna['SE1_FilteredSeqCount'] + \
    metadata_rna['SE2_FilteredSeqCount']
metadata_rna['Combined_Filtered_FragmentCount'] = \
    metadata_rna['PE_FilteredSeqCount'] + \
    metadata_rna['SE1_FilteredSeqCount'] + \
    metadata_rna['SE2_FilteredSeqCount']

# final assigned reads/fragments -----------------------------------------------

metadata_rna['Combined_AssignedReads_Exon'] = \
    metadata_rna['Exon_PE_Assigned']*2 + \
    metadata_rna['Exon_SE1_Assigned'] + \
    metadata_rna['Exon_SE2_Assigned']

metadata_rna['Combined_AssignedFragments_Exon'] = \
    metadata_rna['Exon_PE_Assigned'] + \
    metadata_rna['Exon_SE1_Assigned'] + \
    metadata_rna['Exon_SE2_Assigned']

metadata_rna['Combined_AssignedReads_Gene'] = \
    metadata_rna['Gene_PE_Assigned']*2 + \
    metadata_rna['Gene_SE1_Assigned'] + \
    metadata_rna['Gene_SE2_Assigned']

metadata_rna['Combined_AssignedFragments_Gene'] = \
    metadata_rna['Gene_PE_Assigned'] + \
    metadata_rna['Gene_SE1_Assigned'] + \
    metadata_rna['Gene_SE2_Assigned']


# final RNA metadata ------------------—------------------—---------------------

metadata_rna.to_csv("Metadata/A08b_metadata_RNA.tsv", sep = "\t")
