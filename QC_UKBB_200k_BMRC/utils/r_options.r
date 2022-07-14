# 03_2_initial_sample_qc_plot.r
PLOTS <- '/well/lindgren/dpalmer/ukb_exome_qc/plots/'

# Define some thresholds

# 03_2_initial_sample_qc_filter.r
T_sample_callRate <- 0.95
T_dpMean <- 19.5
T_gqMean <- 47.8

# 06_1_impute_sex_superpopulations_plot.r
T_impute_sex_lower <- 0.2
T_impute_sex_upper <- 0.8

# 08_2_final_variant_qc_filter.r
T_variant_callRate <- 0.97
T_pHWE <- 1e-10

# 09_2_final_sample_qc_filter.r
n_mads <- 4*1.4826 # Approx 4 sds if the dist were normal
n_mads_singleton <- 20
