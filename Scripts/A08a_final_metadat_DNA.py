
# A08b_final_metadat_DNA.py ====================================================
# assumes no changes to script output names from A05/A07

# setup ------------------------------------------------------------------------

import pandas as pd
import glob


# load tables ------------------------------------------------------------------

def read_tbl_wrapper(filepath, prefix = ""):
    return(pd.read_csv(filepath, delimiter = "\t", index_col = 0).add_prefix(prefix))

df_trim_fastp = read_tbl_wrapper("Metadata/A05a_trimming.tsv", "Premap_")

df_bismark_pe = read_tbl_wrapper("Metadata/A05b_DNA_maprate_PE.tsv", "PE_")
df_bismark_se1trim = read_tbl_wrapper("Metadata/A05b_DNA_maprate_SE1trim.tsv", "SE1t_")
df_bismark_se2trim = read_tbl_wrapper("Metadata/A05b_DNA_maprate_SE2trim.tsv", "SE2t_")
df_bismark_se1unmap = read_tbl_wrapper("Metadata/A05b_DNA_maprate_SE1unmap.tsv", "SE1u_")
df_bismark_se2unmap = read_tbl_wrapper("Metadata/A05b_DNA_maprate_SE2unmap.tsv", "SE2u_")

df_picard_pe = read_tbl_wrapper("Metadata/A05c_DNA_picard_PE.tsv", "PE_")
df_picard_se = read_tbl_wrapper("Metadata/A05c_DNA_picard_SE.tsv", "SE_")

df_DNAfrac = read_tbl_wrapper("Metadata/A05d_DNA_global_mCfracs.tsv", "")

df_samstat_pe = read_tbl_wrapper("Metadata/A05e_DNA_samstats_PE.tsv", "PE_")
df_samstat_se = read_tbl_wrapper("Metadata/A05e_DNA_samstats_SE.tsv", "SE_")

df_coverage_1x = read_tbl_wrapper("Metadata/A05f_DNA_cov_chrXdivY.tsv", "")
df_coverage_sex = read_tbl_wrapper("Metadata/A05f_DNA_cov_percent1x.tsv", "")


# merge on rowname index (wellprefix) ------------------------------------------

# SE1trim, SE2trim = singletons from trimming
# SE1unmap, SE2unmap = singletons from PE-mapping 

metadata_DNA = pd.concat([df_trim_fastp,
                         df_bismark_pe, df_bismark_se1trim, df_bismark_se2trim, df_bismark_se1unmap, df_bismark_se2unmap,
                         df_picard_pe, df_picard_se,
                         df_DNAfrac,
                         df_samstat_pe, df_samstat_se,
                         df_coverage_1x, df_coverage_sex], axis = 1)
metadata_DNA = metadata_DNA.apply(pd.to_numeric, errors = "coerce")



# few 'combined' PE & SE metadata stats ----------------------------------------
# some attempts at combined/weighted metrics

# notes: (i) most output above reports paired-end read-pairs/fragments versus reads, hence *2
#        (ii) reads that are attempted to be remapped in "stage-two" SE mapping
#            are subsets of "PE_TotalReadsIn"; thus can't simply sum(total mapped) / sum(reads in)

# estimated total mapping rate (at read level) ----------------------------------------

metadata_DNA['Combined_TotalReadsIn'] = \
    metadata_DNA['PE_TotalReadPairsIn']*2 + \
    metadata_DNA['SE1t_TotalReadsIn'] + \
    metadata_DNA['SE2t_TotalReadsIn']

metadata_DNA['Combined_UniqueMappedReads'] = \
    metadata_DNA.filter(regex = "UniqueMapped", axis = 1
    ).assign(PE_UniqueMappedPairs=metadata_DNA['PE_UniqueMappedPairs']*2
    ).sum(axis = 1)

metadata_DNA['Combined_ReadMappingRate'] = \
    metadata_DNA['Combined_UniqueMappedReads'] / metadata_DNA['Combined_TotalReadsIn']


# percent of dedupe reads that pass filtering ----------------------------------------

metadata_DNA['PE_Filtered_Rate'] = \
    metadata_DNA['PE_FilteredSeqCount'] / \
    (2*metadata_DNA['PE_picard_npairsin']*(1 - metadata_DNA['PE_picard_perc_dupe']))
metadata_DNA['SE_Filtered_Rate'] = \
    metadata_DNA['SE_FilteredSeqCount'] / \
    (metadata_DNA['SE_picard_nreadsin']*(1 - metadata_DNA['SE_picard_perc_dupe']))
metadata_DNA['Combined_Filtered_ReadCount'] = \
    metadata_DNA['PE_FilteredSeqCount'] + metadata_DNA['SE_FilteredSeqCount']
metadata_DNA['Combined_Filtered_ReadPercent'] = \
    metadata_DNA['Combined_Filtered_ReadCount'] / \
    (metadata_DNA['PE_picard_npairsin']*(1 - metadata_DNA['PE_picard_perc_dupe'])*2 +
    metadata_DNA['SE_picard_nreadsin'])


# final total # fragments --------------------------------------------------------------
# (whereas Combined_Filtered_ReadCount above in terms of reads)

metadata_DNA['Combined_Filtered_FragmentCount'] = \
    metadata_DNA['PE_FilteredSeqCount']/2 + metadata_DNA['SE_FilteredSeqCount']


# final DNA metadata -----------------------------------------------------------

metadata_DNA.to_csv("Metadata/A08a_metadata_DNA.tsv", sep = "\t")
