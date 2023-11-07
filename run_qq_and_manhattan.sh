#!/bin/sh

# Load R module
module load R

phenotypes=("Age_related_macular_degeneration")
pops=("AFR" "AMR" "EAS" "EUR" "SAS")

for phenotype in "${phenotypes[@]}"
do
   for pop in "${pops[@]}"
   do
		sbatch --export=pheno=bam qq_and_manhattan_template.sh
   done
done
