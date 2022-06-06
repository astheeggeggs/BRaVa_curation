library(dplyr)
library(data.table)

source("utils/r_options.r")
source("utils/pretty_plotting.r")

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("--combined_variant_qc_file", required=TRUE, help="Path to COMBINED_VARIANT_QC_FILE")
parser$add_argument("--variant_qc_list", required=TRUE, help="Path to VARIANT_LIST prefix to output for final variant filtering")
args <- parser$parse_args()

COMBINED_VARIANT_QC_FILE <- args$combined_variant_qc_file
VARIANT_QC_LIST <- args$variant_qc_list
VARIANT_SUMMARY <- args$variant_summary 

# COMBINED_VARIANT_QC_FILE <- "/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/variants/08_final_qc.variants_combined.tsv"
# VARIANT_QC_LIST <- "/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/variants/08_final_qc.pop.keep.variant_list"
# VARIANT_SUMMARY <- "/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/variants/08_variant_count.pop.tsv"

dt <- fread(COMBINED_VARIANT_QC_FILE, header=TRUE, sep='\t')
print("Initial number of variants:")
print(dt %>% group_by(pop) %>% summarise(n()))

dt_out <- dt %>% filter((variant_qc.AF > 0) & (variant_qc.AF < 1))

# Print, grouped by 1000 genomes label.
print("After removing invariant sites in the cleaned dataset:")
print(dt_out %>% group_by(pop) %>% summarise(n()))

dt_out <- dt_out %>% filter(variant_qc.call_rate >= T_variant_callRate)
print("After ensuring that the overall call rate is good:")
print(dt_out %>% group_by(pop) %>% summarise(n()))

dt_out <- dt_out %>% filter(variant_qc.p_value_hwe > T_pHWE)
print("After ensuring that sites pass HWE p-value threshold:")
print(dt_out %>% group_by(pop) %>% summarise(n()))

dt_final_variant_summary <- data.table(
	rbind(
		dt %>% 
			group_by(Population=pop) %>% summarise(Variants=n()) %>% 
			mutate(filter="Variants after initial filter"),
		dt %>% filter((variant_qc.AF <= 0) | (variant_qc.AF >= 1)) %>% 
			group_by(Population=pop) %>% summarise(Variants=n()) %>% 
			mutate(filter="Invariant sites after sample filters"),
		dt %>% filter(variant_qc.AF > 0, variant_qc.AF < 1, variant_qc.call_rate < T_variant_callRate) %>% 
			group_by(Population=pop) %>% summarise(Variants=n()) %>% 
			mutate(filter=paste0("Overall variant call rate < ", T_variant_callRate)),
		dt %>% filter(variant_qc.AF > 0, variant_qc.AF < 1, variant_qc.p_value_hwe <= T_pHWE) %>% 
			group_by(Population=pop) %>% summarise(Variants=n()) %>% 
			mutate(filter="Variants failing HWE filter"),
		dt_out %>% group_by(Population=pop) %>% summarise(Variants=n()) %>%	
			mutate(filter="Variants after filters")
		)
	)

fwrite(dt_final_variant_summary, file=VARIANT_SUMMARY, quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')
fwrite(dt_out %>% select("pop", "locus", "alleles"), file=VARIANT_QC_LIST, quote=FALSE, row.names=FALSE, col.names=TRUE, sep='\t')
