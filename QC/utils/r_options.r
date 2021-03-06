# Plotting directory for all figures
PLOTS <- '/set/your/plotting/directory'

# Define some thresholds
# 03_2_initial_sample_qc_filter.r
T_sample_callRate <- 0.93
T_pct_contamination <- 0.02
T_pct_chimeras <- 0.015
T_dpMean <- 30
T_gqMean <- 55

# 06_1_impute_sex_plot.r
T_impute_sex_lower <- 0.2
T_impute_sex_upper <- 0.8

# 08_final_variant_qc_filter.r
T_variant_callRate <- 0.97
T_pHWE <- 1e-10