
# ==============================================================================
# Scripts/A01c_well_filepaths.py
# expands plate-level metadata (A01b) into well-level metadata
# ==============================================================================

# recommend running interactively in python/Jupyter to check outputs,
# but shouldn't require any changes to defaults

# load packages ----------------------------------------------------------------

import itertools
import pandas as pd
import numpy as np
import os

# if running interactively, check snmCT_parameters.env loaded or manually spec os.environ e.g.,
# os.environ['projdir'] ="/u/project/cluo/chliu/Analyses/IGVF"; os.chdir(os.environ['projdir'])
# os.environ['metadat_plate'] = "Metadata/A01b_plate_metadata.csv"
# os.environ['metadat_well'] = "Metadata/A01c_well_filepath.csv"



# expand A01b metadata by well -------------------------------------------------

# load A01b
plates_df = pd.read_csv(os.environ['metadat_plate'], index_col=0)

# from pandas documentation
def expand_grid(data_dict):
    """Create a dataframe from every combination of given values."""
    rows = itertools.product(*data_dict.values())
    return pd.DataFrame.from_records(rows, columns=data_dict.keys())

filepath_df = expand_grid({'plate': plates_df['plate'],
    'row' : [chr(x) for x in range(65, 65+16)],
    'col' : [str(x + 1) for x in range(24)]})
filepath_df['well'] = filepath_df[['row', 'col']].agg(''.join, axis = 1)
filepath_df['wellprefix'] = filepath_df['plate'] + "_" + filepath_df['well']

filepath_df = pd.merge(filepath_df, plates_df, how = "left", on = "plate")



# batch into sets of 24 for bismark, STAR processing steps ---------------------
# (by default, one row at a time, incremented by platenum)

# - alternatively, could make smaller batches of wells (e.g., n = 5) for compute
#   environments that favor many small jobs versus a few long jobs,
# - or two sets of batches e.g., filepath_df['batchnum_A04a_bismark']
#   pulled by the sub scripts for the A04a script only

nwellstot = filepath_df.shape[0]
wells_per_batch = 24
filepath_df['batchnum'] =\
    pd.Series(range(0, np.ceil(nwellstot / wells_per_batch).astype(int))
             ).repeat(wells_per_batch)[0:nwellstot].reset_index(drop = True) + 1

print( "number of total wells:" )
print( nwellstot )

filepath_df.index = filepath_df.index.astype(int) + 1

def basename(pathin):
    return(pathin.split("/")[-1])

print( "number of plates:" )
print( "Nplates: " + str( filepath_df['platenum'].max() ) )

print( "number of batches:" )
print( "Nbatches: " + str( filepath_df['batchnum'].max() ) )



# then extensive file paths for sections A02-A06 -------------------------------
# (inelegant, but useful for file checking/compiling info)

# A02: demultiplexing 
# all in dir: fastq_demultip/

filepath_df['A02a_fqgz_demultip_R1'] = "fastq_demultip/" + filepath_df[['plate', 'well']].agg('_'.join, axis = 1) + "_indexed_R1.fastq.gz"
filepath_df['A02a_fqgz_demultip_R2'] = "fastq_demultip/" + filepath_df[['plate', 'well']].agg('_'.join, axis = 1) + "_indexed_R2.fastq.gz"

filepath_df['A02a_txt_summary1'] = "fastq_demultip/" + filepath_df['plate'] + "_summary_1.txt"
filepath_df['A02a_txt_summary2'] = "fastq_demultip/" + filepath_df['plate'] + "_summary_2.txt"



# A03: trimming ----------------------------------------------------------------
# all in dir: fastq_trimmed/

filepath_df['A03a_fqgz_paired_R1'] = "fastq_trimmed/" + filepath_df['wellprefix'] + "_paired_R1.fastq.gz"
filepath_df['A03a_fqgz_paired_R2'] = "fastq_trimmed/" + filepath_df['wellprefix'] + "_paired_R2.fastq.gz"

filepath_df['A03a_fqgz_singletrim_R1'] = "fastq_trimmed/" + filepath_df['wellprefix'] + "_singletrim_R1.fastq.gz"
filepath_df['A03a_fqgz_singletrim_R2'] = "fastq_trimmed/" + filepath_df['wellprefix'] + "_singletrim_R2.fastq.gz"

filepath_df['A03a_json_fastp'] = "fastq_trimmed/" + filepath_df['wellprefix'] + ".json"



# A04: bismark -----------------------------------------------------------------

filepath_df['A04a_dir_bismark'] = "mapping_bismark/" + filepath_df['wellprefix'] + "/"

# (i) paired-end mapping outputs
filepath_df['A04a_bam_bismark_PE'] = \
filepath_df['A04a_dir_bismark'] + filepath_df['A03a_fqgz_paired_R1'].apply(basename).str.replace(".fastq.gz", "_bismark_bt2_pe.bam")
filepath_df['A04a_fqgz_unmap_R1'] = \
filepath_df['A04a_dir_bismark'] + filepath_df['A03a_fqgz_paired_R1'].apply(basename) + "_unmapped_reads_1.fq.gz"
filepath_df['A04a_fqgz_unmap_R2'] = \
filepath_df['A04a_dir_bismark'] + filepath_df['A03a_fqgz_paired_R2'].apply(basename) + "_unmapped_reads_2.fq.gz"

# single-end mapping outputs
filepath_df['A04a_bam_bismark_SE1trim'] = filepath_df['A04a_dir_bismark'] + filepath_df['A03a_fqgz_singletrim_R1'].apply(basename).str.replace(".fastq.gz", "_bismark_bt2.bam")
filepath_df['A04a_bam_bismark_SE2trim'] = filepath_df['A04a_dir_bismark'] + filepath_df['A03a_fqgz_singletrim_R2'].apply(basename).str.replace(".fastq.gz", "_bismark_bt2.bam")

filepath_df['A04a_bam_bismark_SE1unmap'] = filepath_df['A04a_dir_bismark'] + filepath_df['A04a_fqgz_unmap_R1'].str.replace(".fq.gz", "_bismark_bt2.bam")
filepath_df['A04a_bam_bismark_SE2unmap'] = filepath_df['A04a_dir_bismark'] + filepath_df['A04a_fqgz_unmap_R2'].str.replace(".fq.gz", "_bismark_bt2.bam")

# bismark logs
filepath_df['A04a_txt_bismark_PE'] = filepath_df['A04a_dir_bismark'] +\
filepath_df['wellprefix'] + "_paired_R1_bismark_bt2_PE_report.txt"
filepath_df['A04a_txt_bismark_SE1unmap'] = filepath_df['A04a_dir_bismark'] +\
filepath_df['wellprefix'] + "_paired_R1.fastq.gz_unmapped_reads_1_bismark_bt2_SE_report.txt"
filepath_df['A04a_txt_bismark_SE2unmap'] = filepath_df['A04a_dir_bismark'] +\
filepath_df['wellprefix'] + "_paired_R2.fastq.gz_unmapped_reads_2_bismark_bt2_SE_report.txt"
filepath_df['A04a_txt_bismark_SE1trim'] = filepath_df['A04a_dir_bismark'] +\
filepath_df['wellprefix'] + "_singletrim_R1_bismark_bt2_SE_report.txt"
filepath_df['A04a_txt_bismark_SE2trim'] = filepath_df['A04a_dir_bismark'] +\
filepath_df['wellprefix'] + "_singletrim_R2_bismark_bt2_SE_report.txt"

# (ii) picard de-duplication
filepath_df['A04a_bam_dedupe_PE'] = filepath_df['A04a_dir_bismark'] + "PE_dedupe.bam"
filepath_df['A04a_bam_merge_SE'] = filepath_df['A04a_dir_bismark'] + "SE_merge.bam"
filepath_df['A04a_bam_mergesort_SE'] = filepath_df['A04a_dir_bismark'] + "SE_mergesort.bam"
filepath_df['A04a_bam_mergesortdedupe_SE'] = filepath_df['A04a_dir_bismark'] + "SE_mergesortdedupe.bam"

filepath_df['A04a_txt_picard_PE'] = filepath_df['A04a_dir_bismark'] + "picard_PE.log"
filepath_df['A04a_txt_picard_SE'] = filepath_df['A04a_dir_bismark'] + "picard_SE.log"

# (iii) read-level filtering
filepath_df['A04a_sam_dedupeq10_PE'] = filepath_df['A04a_dir_bismark'] + "PE.dedupe_q10.sam"
filepath_df['A04a_sam_dedupeq10_SE'] = filepath_df['A04a_dir_bismark'] + "SE.dedupe_q10.sam"

filepath_df['A04a_allc_final'] = filepath_df['A04a_dir_bismark'] + "allc.tsv.gz"
filepath_df['A04a_allctbi_final'] = filepath_df['A04a_dir_bismark'] + "allc.tsv.gz.tbi"
filepath_df['A04a_txt_allccheck'] = filepath_df['A04a_dir_bismark'] + "allc_check.txt"

# sam stats for coverage, final counts
filepath_df['A04e_txt_samstats_PE'] = filepath_df['A04a_dir_bismark'] + "samstats_PE"
filepath_df['A04e_txt_samstats_SE'] = filepath_df['A04a_dir_bismark'] + "samstats_SE"

filepath_df['A04f_txt_covtot'] = filepath_df['A04a_dir_bismark'] + "total_cov_by_chr"
filepath_df['A04f_txt_covnsites'] = filepath_df['A04a_dir_bismark'] + "nbases_cov_by_chr"



# A05: STAR mapping ------------------------------------------------------------

filepath_df['A05a_dir_star'] = "mapping_star/" + filepath_df['wellprefix'] + "/"

# paired-end mapping outputs (A05a)
filepath_df['A05a_bam_star_PE'] = filepath_df['A05a_dir_star'] + "PE.Aligned.out.bam"
filepath_df['A05a_bam_star_SE1'] = filepath_df['A05a_dir_star'] + "SE1.Aligned.out.bam"
filepath_df['A05a_bam_star_SE2'] = filepath_df['A05a_dir_star'] + "SE2.Aligned.out.bam"

filepath_df['A05a_fq_unmap_R1'] = filepath_df['A05a_dir_star'] + "PE.Unmapped.out.mate1"
filepath_df['A05a_fq_unmap_R2'] = filepath_df['A05a_dir_star'] + "PE.Unmapped.out.mate2"

filepath_df['A05a_txt_star_PE'] = filepath_df['A05a_dir_star'] + "PE.Log.final.out"
filepath_df['A05a_txt_star_SE1'] = filepath_df['A05a_dir_star'] + "SE1.Log.final.out"
filepath_df['A05a_txt_star_SE2'] = filepath_df['A05a_dir_star'] + "SE2.Log.final.out"

# filtered outputs (A05c)
filepath_df['A05c_bam_starfilt_PE'] = filepath_df['A05a_dir_star'] + "PE.Final.bam"
filepath_df['A05c_bam_starfilt_SE1'] = filepath_df['A05a_dir_star'] + "SE1.Final.bam"
filepath_df['A05c_bam_starfilt_SE2'] = filepath_df['A05a_dir_star'] + "SE2.Final.bam"

# samtools & picard output (A05e)
filepath_df['A05e_txt_samtools_PE'] = filepath_df['A05a_dir_star'] + "samstats_PE"
filepath_df['A05e_txt_samtools_SE1'] = filepath_df['A05a_dir_star'] + "samstats_SE1"
filepath_df['A05e_txt_samtools_SE2'] = filepath_df['A05a_dir_star'] + "samstats_SE2"

filepath_df['A05e_txt_picard_PE'] = filepath_df['A05a_dir_star'] + "picard_PE"
filepath_df['A05e_txt_picard_SE1'] = filepath_df['A05a_dir_star'] + "picard_SE1"
filepath_df['A05e_txt_picard_SE2'] = filepath_df['A05a_dir_star'] + "picard_SE2"



# finally, export --------------------------------------------------------------

print(filepath_df.shape)
filepath_df.to_csv(os.environ['metadat_well'])



