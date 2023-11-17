#!/bin/bash

#SBATCH -J annot_long 
#SBATCH -o annot_long-%j.out 
#SBATCH -e annot_long-%j.err 
#SBATCH -p short 

echo "------------------------------------------------" 
echo "Run on host: "`hostname` 
echo "Operating system: "`uname -s` 
echo "Username: "`whoami` 
echo "Started at: "`date` 
echo "------------------------------------------------" 

# Load R module
module purge
module load R

spliceai_dir="/well/lindgren/barney/spliceai/out"
vep_annotations_dir="/well/lindgren/barney/variant-annotation/run-vep/out"
out_dir="/well/lindgren/UKBIOBANK/dpalmer/annotation_munging"
annotation_repo_dir="/well/lindgren/UKBIOBANK/dpalmer/variant-annotation"
mkdir -p ${out_dir}

for chr in {{1..22},"X"}
do
    vep="${vep_annotations_dir}/ukb_wes_450k.qced.chr${chr}.vep_processed.txt"
    spliceai="${spliceai_dir}/ukb_wes_450k.qced.v6.sites_only.${chr}.all.vcf"
    out="${out_dir}/ukb_wes_450k.july.qced.brava_common_rare.v7.chr${chr}.saige.txt"
    python3 ${annotation_repo_dir}/brava_create_annot.py --vep ${vep} --spliceai ${spliceai} --out_file $out
done
