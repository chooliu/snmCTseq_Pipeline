# Public Release, v1


# version 1.2 

## v1.2.1

Minor fix: commenting out parts of A00d, A01b, A01c that should only be run if executing interactively (code chunk starting with "if running interactively [...] os.environ [...]")


## v1.2.0

**Major changes:**
- Added `Documentation/submission_helper.txt` which collects all the commands to run in one place. Namely, all the submission scripts (.sub) with helper `-t 1-Nplates` and `-t 1-Nbatches` info.
- Demultiplexing by default now incorporates 1bp fuzzy matching for cell barcode sequences. Because each cell barcode is at least 2bp different from every other cell barcode (Hamming distance >= 2), by default a 1bp mismatch is allowed. This tends to recover ~0.5 to 1.0% more read pairs in exchange for some extra computational time.  This behavior can be disabled by changing` $enable_fuzzy=1` to `=0` in `A02a_demultiplex_fastq.pl`. 
- Rearrange order of notebooks to group modalities. i.e., A04+A05 = methylation and A06+A07 = transcriptome, as opposed to previous A04+A05=mapping and A06+A07=metadata generation, (alternating mC and RNA). Changes many filepaths. 
- Split former comprehensive "`A04a_bismark_map_to_allc.sub`" script into three parts (A04a [map], A04b [filter reads], A04c [allc]) for far easier cluster scheduling, albeit at the expense of more complex job dependencies throughout the A04* series. The all-in-one previously requested a 48-hour queue that didn't play well with available UCLA queues.
- Also split RNA-mapping via start into a separate mapping and filtering step.
- Changed default behavior of RNA deduplication from `--DUPLICATE_SCORING_STRATEGY` from random to instead keeping the highest quality read (`Picard MarkDuplicates` default).

**Minor changes:**

- A00/A01 (setup/metadata):
    - change example to lambda genome since we now do lambda spike-in by default for mCT
    - change example study displayed (32-plate IGVF sequencing run; consortium-specified hg38 used)
    - "`envvar_needed`"/"`envspec`" variables included to extract from `snmCTseq.env` to limit manual copying in interactive scripts like A01d, A02a, A05, A07
- A02+A03 (demultip/trimming): 
    - FastQC bugs in last version due to some variables  `$numwells_run`.
    - clean up all small scripts e.g., "A03b_check_trimmed" that loops through anticipated output files and checks if the files are available. 
- A04+A05 (methylation):
    - remove "allc_check" file.
    - mcds generation requires a .tsv of .allc names and file locations; this was moved to the mcds generation step. 
    - add helper note on 5kb quantifications for mcds (based on the allcools 5kb-pipeline).
    - omit `samtools stats` step.
- A05, A06, A08:
    - metadata aggregation functions were not all robust to missing inputs. now checks if the file exists and skips files as needed.
    - summarizes numbers of missing values and missing filepaths to help investigate whether particular batches were incomplete.

# version 1.1

## v1.1.3

Only substantial change to pipeline: addition of `picard MarkDuplicates` in RNA filtering steps. Note that `--DUPLICATE_SCORING_STRATEGY RANDOM` is for reducing bias in allele-specific mapping applications. (Have been on the fence about using this as default behavior or not).

Miscellaneous bug fixes (mainly still from shift from hardcoded paths towards use of `.env` file). Special thanks to Nasser for help troubleshooting.
- some scripts incorrectly referenced `$projdir` instead of `$dir_proj`.
- In A01a, Picard CreateSequenceDictionary should use .fa as input.
- In A01b, allow up to 8 lanes (L001-L008) instead of 4 for newer Illumina sequencers.
- In A01c, set `regex = True` to future proof `pd.str.replace`. (Default behavior will be False in figure)
- In A01d genome reference setup, `refFlat.txt.gz` generation had columns in incorrect order, resulting in errors downstream using Picard CollectRnaSeqMetrics. (needed to paste column 12, then 1-10; old version had other way around). 
- function in scripts like A04b checking the number of missing files per batch weren't working properly, needed to explicitly re-code input as bash array. to delineate the fixed version, corrected function renamed from `check_filepath_in_batch` into  `check_filepath_by_batch`.
- A06: fix one instance of hard-coded `ref_chromsizes`.
- reverting one A07 change from v1.1.2 where it should have been relative to `metadata_plate`. 
- increase some memory/time requirements
- updated readme, added more FAQs to detailed_overview


## v1.1.2

- Change paths in `A00*` scripts to be relative to `.env` parameters (versus hardcoded `.fasta` names).
- Minor fix to `A04d`: blank lines causing indexing problems when compiling global methylation levels into `Metadata/A04d_mCfrac_*.tsv` (`wellprefix` row-indices missing and looked duplicated to pandas during metadata compile steps). Also header not being written properly when run in "overwrite" mode, which could cause some entries to be read as headers (apparent missing rows/failed vertical join in compile steps).
- Major issue in metadata compile scripts A06*, A07*, should reference `metadata_well` for filepaths not `metadata_plate` (e.g., couldn't load expected mapping rate files).

## v1.1.1

A few critical bugs related to variable names from new v1.1.0 `snmCT_parameters.env` file structure:
- in A01, `$dir_originalfastq` didn't match .env param name
- annotations in ref file should be .bed instead of .tsv (e.g., should be `annotations/genebody.bed` not `genebody.tsv`)
- forgot to add the FastQC randomly selected wells (`$wells_to_run`) of A02c, A03c in `.env` file. 

Minor fixes: 
- explicitly add `$ref_dir/annotations/` for A00* outputs.
- minor qsub jobs name/dependency fixes
- minor formatting updates, added more comments in Notebooks & D
- unfortunately, notebook metadata cleaning not working (e.g., cellid, ipykernel names updated -- ignore these in the commit)

# version 1.0

## v1.1.0
- Added this documentation, public Github in prep for production-scale consortium work.
- Parameters now kept in one environmental variable file (`snmCT_parameters.env`) in lieu of specifying parameters like `$projdir` that were repeated at the start of almost every script. This requires an `xargs --> export` command & `import os` to read environmental variables for certain python scripts.
- Converting some previously interactively run scripts to also be `qsub` submissions for consistency.
- Explicitly put cell barcode `.fa` safelists in Notebook and Scripts for modularity, although the safelists are the same as in [@luogenomics/demultiplexing](https://github.com/luogenomics/demultiplexing).
- Change default of 48 wells/batch &rarr; 24 wells/batch (one row of wellplate) since some difficulty reporting getting larger resource nodes. Scaled time requests accordingly.

## v1.0
- Large structural/organizational changes without other alterations from v0 prototypes: Since the intermediates and final output for each well are fixed, made extensive changes to workflow structure to be based on pre-specified, known output filepaths in a "targets" file, and looping through these filepaths. 
- Some code now clunkier (reliance on `${basharrays[@]}`), but helps with failed job checking and decreases the degree of user input required at each stage.

# Prototypes, v0 

* **v0.5:** Changes to .bam inputs for STAR mCH/CH (although not the underlying filtering script).
    - Some issues with folks using different versions of `subread`/`featureCounts` (tested on v2.0.1). Later versions more stringent, e.g., throw errors instead of warnings if a mixture of paired-end and single-end reads are input for quantification. 
    - Relevant because "mapping singletons" (e.g., R2 aligns but R1 does not) were previously included in the PE-alignments and fed into featureCounts along with proper pairs; unclear based on developer comments if these were correctly quantified. To address this, we explicitly now subset to only proper pairs before filtering and `featureCounts`.
    - Might be other residual "future proofing" that needs to be done for newer featureCounts: e.g., noticed more stringent error checking (failure on empty .bam file in). 
* **v0.4:** Because mC quantification in `allcools` didn't know how to work with paired-end alignments (e.g., strand-flipped positions with guanine in the `.allc` file), early mC mapping scripts created an intermediate `.bam` with edited flags before `bam_to_allc` conversion to make sure the correct strand and sequence context was extracted. The bam flag editing to artificially strand-swap were removed in this version.
* **v0.1-0.3:** Proof-of-concept paired-end alignments and development of "two-stage" mapping approach to deal with potential chimeras. Major checks included whether mC and RNA quantifies properly (e.g., do overlapping reads correctly get assigned a coverage of 1 or of 2?), treatment of chimeric sequences, and modifying mCH/CH filtering steps to take directionalities into account. Benchmarking of _in silico_ separation of DNA and RNA library against DNA-only snmC-seq2 and RNA-only smart-seq2/3 datasets.
