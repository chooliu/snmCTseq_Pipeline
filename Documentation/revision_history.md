
# Revision History
Key updates from past iterations (& some superfluous details as notes for myself)

# Public Release, v1

**v1.1.1**

A few critical bugs related to variable names from new v1.1.0 `snmCT_parameters.env` file structure:
- in A01, `$dir_originalfastq` didn't match .env param name
- annotations in ref file should be .bed instead of .tsv (e.g., should be `annotations/genebody.bed` not `genebody.tsv`)
- forgot to add the FastQC randomly selected wells (`$wells_to_run`) of A02c, A03c in `.env` file. 

Minor fixes: 
- explicitly add `$ref_dir/annotations/` for A00* outputs.
- minor qsub jobs name/dependency fixes
- minor formatting updates, added more comments in Notebooks & D
- unfortunately, notebook metadata cleaning not working (e.g., cellid, ipykernel names updated -- ignore these in the commit)

**v1.1.0**
- Added this documentation, public Github in prep for production-scale consortium work.
- Parameters now kept in one environmental variable file (`snmCT_parameters.env`) in lieu of specifying parameters like `$projdir` that were repeated at the start of almost every script. This requires an `xargs --> export` command & `import os` to read environmental variables for certain python scripts.
- Converting some previously interactively run scripts to also be `qsub` submissions for consistency.
- Explicitly put cell barcode `.fa` safelists in Notebook and Scripts for modularity, although the safelists are the same as in [@luogenomics/demultiplexing](https://github.com/luogenomics/demultiplexing).
- Change default of 48 wells/batch &rarr; 24 wells/batch (one row of wellplate) since some difficulty reporting getting larger resource nodes. Scaled time requests accordingly.

### v1.0
- Large structural/organizational changes without other alterations from v0 prototypes: Since the intermediates and final output for each well are fixed, made extensive changes to workflow structure to be based on pre-specified, known output filepaths in a "targets" file, and looping through these filepaths. 
- Some code now clunkier (reliance on `${basharrays[@]}`), but helps with failed job checking and decreases the degree of user input required at each stage.

# Prototypes, v0 
* **v0.5:** Changes to .bam inputs for STAR mCH/CH (although not the underlying filtering script).
    - Some issues with folks using different versions of `subread`/`featureCounts` (tested on v2.0.1). Later versions more stringent, e.g., throw errors instead of warnings if a mixture of paired-end and single-end reads are input for quantification. 
    - Relevant because "mapping singletons" (e.g., R2 aligns but R1 does not) were previously included in the PE-alignments and fed into featureCounts along with proper pairs; unclear based on developer comments if these were correctly quantified. To address this, we explicitly now subset to only proper pairs before filtering and `featureCounts`.
    - Might be other residual "future proofing" that needs to be done for newer featureCounts: e.g., noticed more stringent error checking (failure on empty .bam file in). 
* **v0.4:** Because mC quantification in `allcools` didn't know how to work with paired-end alignments (e.g., strand-flipped positions with guanine in the `.allc` file), early mC mapping scripts created an intermediate `.bam` with edited flags before `bam_to_allc` conversion to make sure the correct strand and sequence context was extracted. The bam flag editing to artificially strand-swap were removed in this version.
* **v0.1-0.3:** Proof-of-concept paired-end alignments and development of "two-stage" mapping approach to deal with potential chimeras. Major checks included whether mC and RNA quantifies properly (e.g., do overlapping reads correctly get assigned a coverage of 1 or of 2?), treatment of chimeric sequences, and modifying mCH/CH filtering steps to take directionalities into account. Benchmarking of _in silico_ separation of DNA and RNA library against DNA-only snmC-seq2 and RNA-only smart-seq2/3 datasets.
