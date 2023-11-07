#!/bin/sh

# Load R module
module load R

phenotypes=("Age_related_macular_degeneration")
pops=("AFR" "AMR" "EAS" "EUR" "SAS")

for phenotype in "${phenotypes[@]}"
do
   for pop in "${pops[@]}"
   do
		sbatch --export=pheno=bam run_qq_and_manhattan.sh
   done
done
