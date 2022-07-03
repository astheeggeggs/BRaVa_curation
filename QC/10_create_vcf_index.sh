#!/usr/bin/env bash
#$ -N tabix_preparation
#$ -cwd
#$ -o /well/lindgren/dpalmer/logs/
#$ -e /well/lindgren/dpalmer/logs/
#$ -P lindgren.prjc
#$ -q short.qe@@short.hge
#$ -t 1-23

source /well/lindgren/dpalmer/ukb_utils/bash/qsub_utils.sh

CHR=$(get_chr ${SGE_TASK_ID})
TRANCHE='200k'

module purge
module load BCFtools

QC_VCF_PREFIX="/well/lindgren/UKBIOBANK/dpalmer/wes_${TRANCHE}/ukb_wes_qc/data/final_mt/10_strict_filtered_chr${CHR}"

declare -a arr=("AFR" "AMR" "EAS" "EUR" "SAS")

for i in "${arr[@]}"
do
	echo "$i"
	vcfFile="${QC_VCF_PREFIX}.${i}.vcf.bgz"
	tabix -p vcf -C ${vcfFile}
done
