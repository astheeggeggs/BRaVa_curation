library(data.table)
# Plotting functions:
source('utils/pretty_plotting.r')

# Thresholds and plotting file locations defined in r_options_BipEx.r
# These are currently set at values that made sense for UKBB 200k exomes.
# Please these numbers based on what the plotting below looks like!
source("utils/r_options.r")

# This is commented because args are called from 03_2_initial_sample_qc_filter.r, which also runs all this and 
# generates the plots.

# suppressPackageStartupMessages(library("argparse"))
# parser <- ArgumentParser()
# parser$add_argument("--tranche", default='200k', help = "Which exome sequencing tranche?")
# args <- parser$parse_args()

CHR <- 1

# Inputs:
INITIAL_SAMPLE_QC_FILE <- paste0(
    '/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE,
    '/ukb_wes_qc/data/samples/03_chr', CHR, '_initial_sample_qc.tsv.bgz')

# Outputs:
INITIAL_COMBINED_SAMPLE_QC_FILE <- paste0(
    '/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE, '/ukb_wes_qc/data/samples/03_initial_sample_qc.tsv')

parser$add_argument("--initial_sample_qc_file", required=TRUE, help="Path to INITIAL_SAMPLE_QC_FILE output from 03_0_initial_sample_qc.py")
parser$add_argument("--sample_information", required=TRUE,
    help=paste0("Path to sample information file (aka phenotype file) - this should contain two columns PCT_CHIMERAS and ",
        "PCT_CONTAMINATION, the chimeric read % and freemix contamination %, taken from the GATK/picard metadata. Also include any ",
        "factor you would like to split on and edit the commented code in this file to plot.")
args <- parser$parse_args()

INITIAL_SAMPLE_QC_FILE <- args$initial_sample_qc_file
SAMPLE_INFORMATION <- args$sample_information

# Inputs 
# INITIAL_SAMPLE_QC_FILE (output from 03_0_initial_sample_qc.py)
# SAMPLE_INFORMATION (aka the phenotype file - input from the user - what sample information is available?)
# In particular, this file should include the FREEMIX contaminiation and PCT_CHIMERAS (Chimeric read percentaage) information (if available).
# If the GATK pipeline was run, these metrics should be available in the GATK/Picard metadeta.

dt <- fread(cmd = paste('zcat', INITIAL_SAMPLE_QC_FILE),
    stringsAsFactors=FALSE, sep='\t', header=TRUE) %>% select(c(s, starts_with('qc_padded_target'), starts_with('gq'), starts_with('dp')))
dt <- dt %>% mutate(
    dp.stdev_sum = (dp.stdev)^2 * dp.n,
    gq.stdev_sum = (gq.stdev)^2 * gq.n
    )

setkey(dt, 's')

for (CHR in seq(2,22)) {
    # Input files
    cat(paste0("chromosome ", CHR, "\n"))
    INITIAL_SAMPLE_QC_FILE <- paste0(
        '/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE,
        '/ukb_wes_qc/data/samples/03_chr', CHR, '_initial_sample_qc.tsv.bgz')

    dt_tmp <- fread(
        cmd = paste('zcat', INITIAL_SAMPLE_QC_FILE),
        stringsAsFactors=FALSE, sep='\t', header=TRUE) %>% select(c(s, starts_with('qc_padded_target'), starts_with('gq'), starts_with('dp')))
    names(dt_tmp) <- c('s', paste0('tmp_',names(dt_tmp)[-1]))
    setkey(dt_tmp, 's')

    dt  <- merge(dt, dt_tmp)

    dt <- dt %>% mutate(
        # Mins
        qc_padded_target.dp_stats.min = pmin(qc_padded_target.dp_stats.min, tmp_qc_padded_target.dp_stats.min),
        qc_padded_target.gq_stats.min = pmin(qc_padded_target.gq_stats.min, tmp_qc_padded_target.gq_stats.min),
        
        # Maxs
        qc_padded_target.dp_stats.max = pmax(qc_padded_target.dp_stats.max, tmp_qc_padded_target.dp_stats.max),
        qc_padded_target.gq_stats.max = pmax(qc_padded_target.gq_stats.max, tmp_qc_padded_target.gq_stats.max),
        
        # Counts
        qc_padded_target.n_called = qc_padded_target.n_called + tmp_qc_padded_target.n_called,
        qc_padded_target.n_not_called = qc_padded_target.n_not_called + tmp_qc_padded_target.n_not_called,
        qc_padded_target.n_filtered = qc_padded_target.n_filtered + tmp_qc_padded_target.n_filtered,
        qc_padded_target.n_hom_ref = qc_padded_target.n_hom_ref + tmp_qc_padded_target.n_hom_ref,
        qc_padded_target.n_het = qc_padded_target.n_het + tmp_qc_padded_target.n_het,
        qc_padded_target.n_hom_var = qc_padded_target.n_hom_var + tmp_qc_padded_target.n_hom_var,
        qc_padded_target.n_non_ref = qc_padded_target.n_non_ref + tmp_qc_padded_target.n_non_ref,
        qc_padded_target.n_singleton = qc_padded_target.n_singleton + tmp_qc_padded_target.n_singleton,
        qc_padded_target.n_snp = qc_padded_target.n_snp + tmp_qc_padded_target.n_snp,
        qc_padded_target.n_insertion = qc_padded_target.n_insertion + tmp_qc_padded_target.n_insertion,
        qc_padded_target.n_deletion = qc_padded_target.n_deletion + tmp_qc_padded_target.n_deletion,
        qc_padded_target.n_transition = qc_padded_target.n_transition + tmp_qc_padded_target.n_transition,
        qc_padded_target.n_transversion = qc_padded_target.n_transversion + tmp_qc_padded_target.n_transversion,
        qc_padded_target.n_star = qc_padded_target.n_star + tmp_qc_padded_target.n_star,

        # Means
        dp.sum = dp.sum + tmp_dp.sum,
        gq.sum = gq.sum + tmp_gq.sum,

        # Stdevs
        dp.stdev_sum = dp.stdev_sum + (tmp_dp.stdev)^2 * tmp_dp.n,
        gq.stdev_sum = dp.stdev_sum + (tmp_gq.stdev)^2 * tmp_gq.n,

        # N
        dp.n = dp.n + tmp_dp.n,
        gq.n = gq.n + tmp_gq.n
    ) %>% select(c(s, starts_with('qc_padded_target'), starts_with('dp'), starts_with('gq')))
    setkey(dt, 's')
}

dt <- dt %>% mutate(
    # Ratios
    qc_padded_target.call_rate = qc_padded_target.n_called / (qc_padded_target.n_called + qc_padded_target.n_not_called + qc_padded_target.n_filtered),
    qc_padded_target.r_ti_tv = ifelse(qc_padded_target.n_transversion == 0, NA, qc_padded_target.n_transition / qc_padded_target.n_transversion),
    qc_padded_target.r_het_hom_var = ifelse(qc_padded_target.n_hom_var == 0, NA, qc_padded_target.n_het / qc_padded_target.n_hom_var),
    qc_padded_target.r_insertion_deletion = ifelse(qc_padded_target.n_deletion == 0, NA, qc_padded_target.n_insertion / qc_padded_target.n_deletion),
    
    # Means
    qc_padded_target.dp_stats.mean = dp.sum / dp.n,
    qc_padded_target.gq_stats.mean = gq.sum / gq.n,

    # Stdevs
    qc_padded_target.dp_stats.stdev = sqrt(dp.stdev_sum / dp.n),
    qc_padded_target.gq_stats.stdev = sqrt(gq.stdev_sum / gq.n)
    ) %>% select(c('s', starts_with('qc_padded_target')))

setkey(dt, 's')

df_pheno <- fread(SAMPLE_INFORMATION) # Here we assume that chimeric reads and pct contamination have column names 'PCT_CONTAMINATION' and 'PCT_CHIMERAS' respectively.
df <- merge(df_pheno, df,  by='s')

names(df) <- gsub("qc_padded_target\\.", "", names(df))

# Shuffle the rows for plotting
df <- df[sample(nrow(df), replace=FALSE),]

save_figures <- TRUE
# PDFs, no splitting.
create_pretty_hist(df, aes(x=call_rate), 'Call Rate', T_sample_callRate,
    title='Call Rate', save_figure=save_figures, file=paste0(PLOTS,'03_callRate_hist'))
create_pretty_hist(df, aes(x=PCT_CONTAMINATION), 'Contamination', T_pct_contamination,
    binwidth=0.0001, xlim=c(0,0.1), title='% Contamination', save_figure=save_figures, file=paste0(PLOTS,'03_contamination_hist'))
create_pretty_hist(df, aes(x=PCT_CHIMERAS), 'Chimeric Reads', T_pct_chimeras,
    binwidth=0.00002, xlim=c(0,0.02), title='% Chimeric Reads', save_figure=save_figures, file=paste0(PLOTS,'03_chimeras_hist'))
create_pretty_hist(df, aes(x=dp_stats.mean), 'Mean Depth', T_dpMean,
    binwidth=2, xlim=c(10, 150), title='Mean Depth', save_figure=save_figures, file=paste0(PLOTS,'03_dpMean_hist'))
create_pretty_hist(df, aes(x=gq_stats.mean), 'Mean Genotype Quality', T_gqMean,
    binwidth=0.5, xlim=c(20, 70), title='Mean Genotype Quality', save_figure=save_figures, file=paste0(PLOTS,'03_gqMean_hist'))

# CDFs, no splitting.
create_pretty_cumulative(df, aes(call_rate), 'Call Rate', T_sample_callRate,
    xlim=c(0.75,1), title='Call Rate', save_figure=save_figures, file=paste0(PLOTS,'03_callRate_cdf'))
create_pretty_cumulative(df, aes(PCT_CONTAMINATION), 'Contamination', T_pct_contamination,
    xlim=c(0,0.1), title='% Contamination', save_figure=save_figures, file=paste0(PLOTS,'03_contamination_cdf'))
create_pretty_cumulative(df, aes(PCT_CHIMERAS), 'Chimeric Reads', T_pct_chimeras,
    xlim=c(0,0.02), title='% Chimeric Reads', save_figure=save_figures, file=paste0(PLOTS,'03_chimeras_cdf'))
create_pretty_cumulative(df, aes(dp_stats.mean), 'Mean Depth', T_dpMean,
    xlim=c(10,150), title='Mean Depth', save_figure=save_figures, file=paste0(PLOTS,'03_dpMean_cdf'))
create_pretty_cumulative(df, aes(gq_stats.mean), 'Mean Genotype Quality', T_gqMean,
    xlim=c(20,70), title='Mean Genotype Quality', save_figure=save_figures, file=paste0(PLOTS,'03_gqMean_cdf'))

# Here are a collection of plots to look at the distribution of metrics after splitting on batch, collection etc.
# You can edit these to plot based on the phenotypic information that you have in your samples.

# # Edit the code below to perform your own plotting, splitting on different factors.

# legend_batch <- FALSE
# legend_collection <- TRUE
# legend_phenotype <- TRUE
# save_figures <- TRUE
# y_label_batch <- 'Batch'
# y_label_batch <- ''
# titles <- c('Call Rate',
#     '% Contamination',
#     '% Chimeric Reads',
#     'Mean Depth (DP)',
#     'Mean Genotype Quality (GQ)')
# alpha <- 0.8

# print(df %>% group_by(LOCATION) %>% summarise(n()))
# print(df %>% group_by(PROJECT_OR_COHORT) %>% summarise(n()))

# # Plot after splitting on location, and colouring by project
# create_pretty_boxplots(df, aes(x=LOCATION, y=call_rate), aes(color=PROJECT_OR_COHORT),
#     T_sample_callRate, x_label='Call Rate', y_label=y_label_batch, key_label='Batch',
#     xlim=quantile(df$call_rate, c(0.01, 0.99)), legend=legend_batch, title=titles[1], save_figure=save_figures,
#     file=paste0(PLOTS,'03_callRate_by_collection'), n_ticks=5, alpha=alpha)
# create_pretty_boxplots(df, aes(x=LOCATION, y=PCT_CONTAMINATION), aes(color=PROJECT_OR_COHORT),
#     T_pct_contamination, x_label='% Contamination', y_label=y_label_batch, key_label='Batch',
#     xlim=quantile(df$PCT_CONTAMINATION, c(0.01, 0.99)), legend=legend_batch, title=titles[2], save_figure=save_figures,
#     file=paste0(PLOTS,'03_contaminiation_by_collection'), n_ticks=5, alpha=alpha)
# create_pretty_boxplots(df, aes(x=LOCATION, y=PCT_CHIMERAS), aes(color=PROJECT_OR_COHORT),
#     T_pct_chimeras, x_label='% Chimeric Reads', y_label=y_label_batch, key_label='Batch',
#     xlim=quantile(df$PCT_CHIMERAS, c(0.01, 0.99)), legend=legend_batch, title=titles[3], save_figure=save_figures,
#     file=paste0(PLOTS,'03_chimeras_by_collection'), n_ticks=5, alpha=alpha)
# create_pretty_boxplots(df, aes(x=LOCATION, y=dp_stats.mean), aes(color=PROJECT_OR_COHORT),
#     T_dpMean, x_label='Mean Depth', y_label=y_label_batch, key_label='Batch',
#     xlim=quantile(df$dp_stats.mean, c(0.01, 0.99)), legend=legend_batch, title=titles[4], save_figure=save_figures,
#     file=paste0(PLOTS,'03_dpMean_by_collection'), alpha=alpha)
# create_pretty_boxplots(df, aes(x=LOCATION, y=gq_stats.mean), aes(color=PROJECT_OR_COHORT),
#     T_gqMean, x_label='Mean Genotype Quality', y_label=y_label_batch, key_label='Batch',
#     xlim=quantile(df$gq_stats.mean, c(0.01, 0.99)), legend=legend_batch, title=titles[5], save_figure=save_figures,
#     file=paste0(PLOTS,'03_gqMean_by_collection'), alpha=alpha)

# # Plot after splitting by cohort and plotting by sampling location.
# create_pretty_boxplots(df, aes(x=PROJECT_OR_COHORT, y=call_rate), aes(color=LOCATION),
#     T_sample_callRate, y_label='Batch', x_label='Call Rate', key_label='Collection',
#     xlim=quantile(df$call_rate, c(0.01, 0.99)), legend=legend_collection, title=titles[1],
#     save_figure=save_figures, file=paste0(PLOTS,'03_callRate_by_batch'), alpha=alpha)
# create_pretty_boxplots(df, aes(x=PROJECT_OR_COHORT, y=PCT_CONTAMINATION), aes(color=LOCATION),
#     T_pct_contamination, y_label='Batch', x_label='% Contamination', key_label='Collection',
#     xlim=quantile(df$PCT_CONTAMINATION, c(0.01, 0.99)), legend=legend_collection, title=titles[2],
#     save_figure=save_figures, file=paste0(PLOTS,'03_contamination_by_batch'), alpha=alpha)
# create_pretty_boxplots(df, aes(x=PROJECT_OR_COHORT, y=PCT_CHIMERAS), aes(color=LOCATION),
#     T_pct_chimeras, y_label='Batch', x_label='% Chimeric Reads', key_label='Collection',
#     xlim=quantile(df$PCT_CHIMERAS, c(0.01, 0.99)), legend=legend_collection, title=titles[3],
#     save_figure=save_figures, file=paste0(PLOTS,'03_chimeras_by_batch'), alpha=alpha)
# create_pretty_boxplots(df, aes(x=PROJECT_OR_COHORT, y=dp_stats.mean), aes(color=LOCATION),
#     T_dpMean, y_label='Batch', x_label='Mean Depth', key_label='Collection',
#     xlim=quantile(df$dp_stats.mean, c(0.01, 0.99)), legend=legend_collection, title=titles[4],
#     save_figure=save_figures, file=paste0(PLOTS,'03_dpMean_by_batch'), alpha=alpha)
# create_pretty_boxplots(df, aes(x=PROJECT_OR_COHORT, y=gq_stats.mean), aes(color=LOCATION),
#     T_gqMean, y_label='Batch', x_label='Mean Genotype Quality', key_label='Collection',
#     xlim=quantile(df$gq_stats.mean, c(0.01, 0.99)), legend=legend_collection, title=titles[5],
#     save_figure=save_figures, file=paste0(PLOTS,'03_gqMean_by_batch'), alpha=alpha)

