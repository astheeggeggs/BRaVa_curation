# Set of scripts to generate annotation summary counts and plots

Before starting, create the conda environment to run the R scripts in steps 2 and 3.

```
conda config --get channels
conda config --add channels 'r'           # lowest priority
conda config --add channels 'defaults'
conda config --add channels 'conda-forge' # highest priority

conda env create --file environment.yml
source activate annot_counts
```

## Step 1: Determine AC for all variants
Example shell script to acheive this using plink 1.9 and `--freqx` is present in `01_get_AC.sh`.
Note that the following scripts require the `.frqx` [output format](https://www.cog-genomics.org/plink/1.9/formats#frqx), so if using bcftools etc, the resultant file will need to be munged to match the column headers of `.frqx`.
(CHR, SNP, A1, A2, C(HOM A1), C(HET), C(HOM A2), C(HAP A1), C(HAP A2), C(MISSING)).

## Step 2a: Generate the BRaVa annotations
The files required are generated in Steps 1 and 2 of the [variant-annotation repo](https://github.com/BRaVa-genetics/variant-annotation). 

First, using files `sites_only_output_chr${chr}_vep.processed.txt` (from [step 1 of variant-annnotiotion repo](https://github.com/BRaVa-genetics/variant-annotation#1-run-vep-version-105-with-loftee-v104_grch38)) and your spliceAI annotated sites only VCF from [step 2 of the variant-annotation repo](https://github.com/BRaVa-genetics/variant-annotation#2-run-spliceai) (Note - both of these files should *not* have had the popmax 1% gnomAD filter applied). 

If you don't have `sites_only_output_chr${chr}_vep.processed.txt` files, run
```
bcftools +split-vep out/sites_only_output_chr${chr}_vep.vcf.gz -d -f '%CHROM:%POS:%REF:%ALT %Gene %LoF %REVEL_score %CADD_phred %Consequence %Feature %MANE_SELECT %CANONICAL %BIOTYPE\n' -o out/sites_only_output_chr${chr}_vep.processed.txt
sed -i '1i SNP_ID GENE LOF REVEL_SCORE CADD_PHRED CSQ TRANSCRIPT MANE_SELECT CANONICAL BIOTYPE' out/sites_only_output_chr${chr}_vep.processed.txt
```
to create them.

Next, run the [BRaVa annotation script](https://github.com/BRaVa-genetics/variant-annotation/blob/main/SAIGE_annotations/scripts/brava_create_annot.py) to extract variant annotations for these files:

```
brava_create_annot.py -v vep_table -s spliceai_vcf -w output_file
rm output_file
```

This python script will create two output files, one for input into SAIGE (which we delete as this set retains variants with popMAX > 1% which we don't want to inadvertantly include in any testing), and another named `${output_file}.long.csv.gz`.
Example shell scripts looping over chromosomes are [here](https://github.com/BRaVa-genetics/variant-annotation/blob/main/SAIGE_annotations/scripts/brava_create_annot.sh).

## Step 2: Merge the AC and annotation information in a summary
This script reads in the output of *Step 1* and *Step 2a* to generate four summary files for each call to the script.

```
Rscript 02_merged_annots_and_AC.r --AC_path ${step_1_output.frqx.gz} --vep_spliceAI_processed ${step_2a_output.long.csv.gz} --spliceAI_bins --out ${out}
```
It's likely that you will generate separate summary files for each chromosome and for each population label. See [`02_get_merged_annots_and_AC.sh`](https://github.com/astheeggeggs/BRaVa_curation/blob/main/QC/annotation_summary/02_get_merged_annots_and_AC.sh) for an example using the DNANexus RAP.
These files should be uploaded to your biobank/cohort specific bucket. We recommend naming the files with the following naming convention so that *Step 3* plotting can be easily carried out:

```
out={biobank}.{pop}.chr{chr}
```

## Step 3 (optional): Plot the AC and variant counts, split by MAF, MAC and spliceAI bin
We've also included a plotting function `03_count_plotting.r` to visualise the data across chromosomes.

Note that this function assumes that the summary files created in *Step 2* named in the following way:

`{biobank}.{pop}.chr{chr}.BRaVa_annotations_{variant, transcript}_summary.tsv.gz"` and placed in a folder (named `{count_dir}`, say).

To generate a `.pdf` of plots, run:

```
Rscript 03_count_plotting.r --count_directory ${count_dir} --biobank ${biobank} --out ${out} --spliceAI_bins
```


