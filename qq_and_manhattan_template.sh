#!/bin/bash

#SBATCH -J plot 
#SBATCH -o plot-%j.out 
#SBATCH -e plot-%j.err 
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

Rscript run_qq_and_manhattan.r --phenotype=$phenotype --population=$pop --sex=$sex
