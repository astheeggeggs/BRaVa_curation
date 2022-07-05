#!/usr/bin/env bash

# UKB genotype calls
ukb_bed="/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_cal_chrX_v2.bed"
ukb_bim="/well/lindgren/UKBIOBANK/DATA/CALLS/ukb_snp_chrX_v2.bim"
ukb_fam="/well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/ukb11867_cal_chr1_v2_s488363.fam"

out="/well/lindgren/UKBIOBANK/dpalmer/ukb_genotype_plink/ukb_snp_chrX_pruned"

# Create the vcf fam
TRANCHE='200k'
pheno="/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/${TRANCHE}/QC_phenotypes.tsv.gz"
vcf_samples="/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/${TRANCHE}/UKBB_WES${TRANCHE}.txt"
zcat $pheno | awk '{print $1, $1}' | tail -n +2 > $vcf_samples

/well/lindgren/dpalmer/plink --bed ${ukb_bed} --bim ${ukb_bim} --fam ${ukb_fam} \
      --keep ${vcf_samples} --maf 0.05 --geno 0.02 \
      --indep 50 5 2 --out ${out}
