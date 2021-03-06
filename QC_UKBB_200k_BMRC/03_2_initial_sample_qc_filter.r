library(dplyr)
library(data.table)
source("utils/r_options.r")

parser <- ArgumentParser()
parser$add_argument("--tranche", default='200k', help = "Which exome sequencing tranche?")
args <- parser$parse_args()

TRANCHE <- args$tranche

# Output
SAMPLE_LIST_INITIAL_QC = paste0('/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE, '/ukb_wes_qc/data/samples/03_initial_qc.keep.sample_list')
SAMPLE_SUMMARY_COUNT = paste0('/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE, '/ukb_wes_qc/data/samples/03_sample_count.tsv')

# Run the plotting again to ensure that the thresholds are as in the plots.
source("03_1_initial_sample_qc_plot.r")

dt_out <- filter(dt, call_rate > T_sample_callRate) %>%
	filter(dp_stats.mean > T_dpMean) %>%
	filter(gq_stats.mean > T_gqMean)

dt_out <- dt_out %>% select(s)
print(dim(dt_out))

fwrite(dt_out, file=SAMPLE_LIST_INITIAL_QC, quote=FALSE, row.names=FALSE, col.names=FALSE)

# Create the table too
dt_summary_count <- data.table(
	"Filter" = c("Samples after initial filter",
			   paste0("Sample call rate < ", T_sample_callRate),
			   paste0("% FREEMIX contamination > ", T_pct_contamination),
			   paste0("% chimeric reads > ", T_pct_chimeras),
			   paste0("Mean DP < ", T_dpMean),
			   paste0("Mean GQ < ", T_gqMean),
			   "Samples after sample QC filters"),
	"Samples" = c(nrow(dt),
			    nrow(filter(dt, call_rate <= T_sample_callRate)),
				nrow(filter(dt, PCT_CONTAMINATION >= T_pct_contamination)),
				nrow(filter(dt, PCT_CHIMERAS >= T_pct_chimeras)),
				nrow(filter(dt, dp_stats.mean <= T_dpMean)),
				nrow(filter(dt, gq_stats.mean <= T_gqMean)),
				nrow(dt_out))
	)

fwrite(dt_summary_count, file=SAMPLE_SUMMARY_COUNT, quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')
