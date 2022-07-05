#!/usr/bin/env bash

# UKB genotype calls
GENO_X_BED=$1
GENO_X_BIM=$2
GENO_X_FAM=$3
PHENOFILE=$4
PRUNED_X_PLINK=$5

# Extract phenotype sample ID from phenotype file to create a fam to restrict to
samples_with_phenotype_data="~/tmp.txt"
zcat $PHENOFILE | awk '{print $1, $1}' | tail -n +2 > $samples_with_phenotype_data

plink --bed ${GENO_BED} --bim ${GENO_BIM} --fam ${GENO_FAM} \
      --keep ${samples_with_phenotype_data} --maf 0.05 --geno 0.02 \
      --indep 50 5 2 --out ${PRUNED_X_PLINK}

rm $samples_with_phenotype_data