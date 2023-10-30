#!/usr/bin/env Rscript

library(argparse)

parser <- ArgumentParser()
parser$add_argument("--phenotype", type="character")
parser$add_argument("--population", type="character")
parser$add_argument("--sex", type="character")
args <- parser$parse_args()

print(paste("Phenotype:", args$phenotype)) 
print(paste("Population:", args$population))
print(paste("Sex:", args$sex))

# create_brava_qq_and_manhattan(save=TRUE, wait_for_completion=FALSE, download_only=TRUE)
create_brava_single_qq_and_manhattan(
	pheno=args$phenotype, pop=args$population, sex=args$sex, save=TRUE
	)
