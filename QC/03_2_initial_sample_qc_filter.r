library(dplyr)
library(data.table)
source("r_functions_and_parameters/r_options_BipEx.r")

# Run the plotting again to ensure that the thresholds are as in the plots.
source("03_initial_sample_qc_plot.r")

df_out <- filter(df, call_rate > T_sample_callRate) %>%
	filter(PCT_CONTAMINATION < T_pct_contamination) %>%
	filter(PCT_CHIMERAS < T_pct_chimeras) %>%
	filter(dp_stats.mean > T_dpMean) %>%
	filter(gq_stats.mean > T_gqMean)

df_out <- df_out %>% select(s)
print(dim(df_out))

fwrite(df_out, file=SAMPLE_LIST_INITIAL_QC, quote=FALSE, row.names=FALSE, col.names=FALSE)

# Create the table too
df_summary_count <- data.table(
	"Filter" = c("Samples after initial filter",
			   paste0("Sample call rate < ", T_sample_callRate),
			   paste0("% FREEMIX contamination > ", T_pct_contamination),
			   paste0("% chimeric reads > ", T_pct_chimeras),
			   paste0("Mean DP < ", T_dpMean),
			   paste0("Mean GQ < ", T_gqMean),
			   "Samples after sample QC filters"),
	"Samples" = c(nrow(df),
			    nrow(filter(df, call_rate <= T_sample_callRate)),
				nrow(filter(df, PCT_CONTAMINATION >= T_pct_contamination)),
				nrow(filter(df, PCT_CHIMERAS >= T_pct_chimeras)),
				nrow(filter(df, dp_stats.mean <= T_dpMean)),
				nrow(filter(df, gq_stats.mean <= T_gqMean)),
				nrow(df_out)),
	"Bipolar Cases" = c(nrow(df %>% filter(PHENOTYPE_COARSE == "Bipolar Disorder")),
			    nrow(filter(df, call_rate <= T_sample_callRate) %>% filter(PHENOTYPE_COARSE == "Bipolar Disorder")),
				nrow(filter(df, PCT_CONTAMINATION >= T_pct_contamination) %>% filter(PHENOTYPE_COARSE == "Bipolar Disorder")),
				nrow(filter(df, PCT_CHIMERAS >= T_pct_chimeras) %>% filter(PHENOTYPE_COARSE == "Bipolar Disorder")),
				nrow(filter(df, dp_stats.mean <= T_dpMean) %>% filter(PHENOTYPE_COARSE == "Bipolar Disorder")),
				nrow(filter(df, gq_stats.mean <= T_gqMean) %>% filter(PHENOTYPE_COARSE == "Bipolar Disorder")),
				length(which(df_out$s %in% (df %>% filter(PHENOTYPE_COARSE == "Bipolar Disorder"))$s))),
	"Controls" = c(nrow(df %>% filter(PHENOTYPE_COARSE == "Control")),
			    nrow(filter(df, call_rate <= T_sample_callRate) %>% filter(PHENOTYPE_COARSE == "Control")),
				nrow(filter(df, PCT_CONTAMINATION >= T_pct_contamination) %>% filter(PHENOTYPE_COARSE == "Control")),
				nrow(filter(df, PCT_CHIMERAS >= T_pct_chimeras) %>% filter(PHENOTYPE_COARSE == "Control")),
				nrow(filter(df, dp_stats.mean <= T_dpMean) %>% filter(PHENOTYPE_COARSE == "Control")),
				nrow(filter(df, gq_stats.mean <= T_gqMean) %>% filter(PHENOTYPE_COARSE == "Control")),
				length(which(df_out$s %in% (df %>% filter(PHENOTYPE_COARSE == "Control"))$s))))

fwrite(df_summary_count, file='../../samples_BipEx/03_sample_count.tsv', quote=FALSE, row.names=FALSE, col.names=FALSE, sep='\t')
