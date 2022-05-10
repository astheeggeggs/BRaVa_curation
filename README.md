# Templates for WES QC pipeline

A collection of scripts to use as a starting point for running quality control on whole exome sequencing datasets using Hail.

Throughout, we assume that the reference geneome is GRCh38. If this isn't the reference used for calling in your dataset, switch to the appropriate reference by changing `REFERENCE = 'GRCh38'` in each of the Hail python scripts.

These scripts assume that you have access to a large number of cores to run them on. If not, you will need to split (e.g. by chromosome) and run scripts in parallel. Examples of this procedure and intermediate plotting (as well as scripts for submitting jobs to a HPC) for very similar steps to these can be found in [this](https://github.com/astheeggeggs/SAIGE_gene_munging/tree/main/QC_scripts) repository.

## Pipeline
### Step 0: Create an initial Hail MatrixTable
Read in joint called VCF and write to Hail MatrixTable format

__Inputs__: 
* Filepath to the jointcalled .vcf to input: `${INPUT_VCF}`

__Outputs__: 
* Filepath to place the output MatrixTable: `${RAW_MT}`

`python 00_load_and_write_vcf_as_mt.py ${INPUT_VCF} ${RAW_MT}`

### Step 1: Genotype QC
* Remove sites with a large number of alleles (>6)
* Filter genotypes based on depth, likelihood, and allele balance.
  * If homozygous reference, at least one of:
    * Genotype quality < 20
    * Depth < 10
  * If heterozygous, at least one of:
    * (Reference allele depth + alternative allele depth)/depth < 0.8
    * (Alternative allele depth)/depth < 0.2
    * Reference phred-scaled genotype posterior < 20
    * Depth < 10
  * If homozygous variant, at least one of:
    * (Alternative allele depth)/depth < 0.8
    * Reference phred-scaled genotype posterior < 20
    * Depth < 10

__Inputs__: 
* Filepath to place the MatrixTable output from step 0: `${RAW_MT}`

__Outputs__:
* Filepath to place the output MatrixTable: `${MT}`
* Filepath to place the output hard-calls MatrixTable: `$MT_HARDCALLS`

`python 01_filterGT.py ${RAW_MT} ${MT} ${MT_HARDCALLS}`

### Step 2: Initial Variant QC

Remove variants that either:
* Fall in a low complexity region
* Fail VQSR (note that if GATK was not used for calling, you will not be able to apply this filter).
* Fall outside of padded target intervals (we recommend 50bp padding)
* Filter out invariant sites

__Inputs__: 
* Filepath to hard-calls MatrixTable from step 1: `${MT_HARDCALLS}`
* Filepath to the target intervals file (see below for expected format): `${TARGET_INTERVALS}`
* Filepath to the padded target intervals file: `${PADDED_TARGET_INTERVALS}`
* Filepath to low-complexity regions intervals file: `${LCRs}`

__Outputs__:
* Filepath to place initial QC metrics for downstream plotting: `${INITIAL_VARIANT_QC_FILE}$
* Filepath to place the set of variants remaining after step 2 filtering: `${INITIAL_VARIANT_LIST}`

`python 02_prefilter_variants.py ${MT_HARDCALLS} ${TARGET_INTERVALS} ${PADDED_TARGET_INTERVALS} ${LCRs} ${INITIAL_VARIANT_QC_FILE} ${INITIAL_VARIANT_LIST}`

We have included an example target intervals file, target intervals file with 50bp padding and LCRs file to show the expected format.

### Step 3: Initial Sample QC

__Inputs__: 
* Filepath to hard-calls MatrixTable from step 1: `${MT_HARDCALLS}`
* Filepath to the set of variants remaining after step 2 filtering: `${INITIAL_VARIANT_LIST}`
* Filepath to the target intervals file (see below for expected format): `${TARGET_INTERVALS}`

__Outputs__:
* Filepath to place initial sample QC metrics for downstream plotting: `${INITIAL_SAMPLE_QC_FILE}$`

`python 03_0_initial_sample_qc.py ${MT_HARDCALLS} ${INITIAL_VARIANT_LIST} ${INITIAL_SAMPLE_QC_FILE}`

### Step 4: High Quality Common Variant Subset
`04_0_export_plink.py`

### Step 5: Determine superpopulation ancestry labels

### Step 6: Sex imputation
This step should be run separately for each superpopulation (MAF differences across superpopulations can throw off the $F$ statistic).
`05_0_impute_sex.py`

# Step 7: Determine related samples
If you have a single homogeneous population, you can use IBD estimation with the following script:
`06_0_ibd.py`

If you have multiple superpopulations, you should use PC-relate to identify related samples with the following script:
`06_0_pc_relate.py`
