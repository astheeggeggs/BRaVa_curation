library(ggplot2)
library(ggsci)
library(data.table)

# Plotting functions:
source('utils/pretty_plotting.r')
# Thresholds and plotting file locations defined in r_options.r
source("utils/r_options.r")

suppressPackageStartupMessages(library("argparse"))

CHR <- 1
TRANCHE <- '200k'
save_figures <- TRUE
PLOTS <- '/well/lindgren/dpalmer/ukb_exome_qc/plots/'

# Input files
SAMPLE_QC_FILE <- paste0("/well/lindgren/UKBIOBANK/dpalmer/wes_", TRANCHE, "/ukb_wes_qc/data/samples/11_sample_metrics_for_plotting_chr", CHR, ".tsv.gz")
VARIANT_QC_FILE <- paste0("/well/lindgren/UKBIOBANK/dpalmer/wes_", TRANCHE, "/ukb_wes_qc/data/variants/11_variant_metrics_for_plotting_chr", CHR, ".tsv.gz")

# Output files
COMBINED_SAMPLE_QC_FILE <- paste0('/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE, '/ukb_wes_qc/data/11_sample_metrics_for_plotting.tsv')
COMBINED_VARIANT_QC_FILE <- paste0('/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE, '/ukb_wes_qc/data/11_variant_metrics_for_plotting.tsv')

dt <- fread(cmd = paste('zcat', VARIANT_QC_FILE), header=TRUE, sep='\t')

dt_list <- list()
dt_list[[1]] <- dt

for (CHR in c(seq(2,22), "X")) {
    # Input files
    cat(paste0("chromosome ", CHR, "\n"))
    VARIANT_QC_FILE <- paste0("/well/lindgren/UKBIOBANK/dpalmer/wes_", TRANCHE, "/ukb_wes_qc/data/variants/11_variant_metrics_for_plotting_chr", CHR, ".tsv.gz")
    dt_tmp <- fread(cmd = paste('zcat', VARIANT_QC_FILE), header=TRUE, sep='\t')
    setkeyv(dt_tmp, c('locus', 'alleles'))
    dt_list[[CHR]] <- dt_tmp
}

dt <- rbindlist(dt_list)

fwrite(dt, file=COMBINED_VARIANT_QC_FILE , sep='\t', quote=FALSE)

dt <- fread(COMBINED_VARIANT_QC_FILE)

# call rate across all variants
create_pretty_hist(dt, aes(x=variant_qc.call_rate), x_label='Call Rate', save_figure=save_figures,
	file=paste0(PLOTS, TRANCHE, '_11_callRate_hist'))
# cumulative call rate
create_pretty_cumulative(dt, aes(x=variant_qc.call_rate), x_label="Call Rate", threshold=NULL,
    key_label='', save_figure=save_figures, file=paste0(PLOTS, TRANCHE, '_11_callRate_cdf'))

# mean DP all variants
create_pretty_hist(dt, aes(x=variant_qc.dp_stats.mean), x_label='Mean DP', save_figure=save_figures,
    file=paste0(PLOTS, TRANCHE, '_11_mean_DP_hist'))
# cumulative mean DP
create_pretty_cumulative(dt, aes(x=variant_qc.dp_stats.mean), x_label="Mean DP", threshold=NULL,
    key_label='', save_figure=save_figures, file=paste0(PLOTS, TRANCHE, '_11_mean_DP_cdf'))

# mean GQ all variants
create_pretty_hist(dt, aes(x=variant_qc.gq_stats.mean), x_label='Mean GQ', save_figure=save_figures,
    file=paste0(PLOTS, TRANCHE, '_11_mean_GQ_hist'))
# cumulative mean GQ
create_pretty_cumulative(dt, aes(x=variant_qc.call_rate), x_label='Mean GQ', threshold=NULL,
    key_label='', save_figure=save_figures, file=paste0(PLOTS, TRANCHE, '_11_mean_GQ_cdf'))

# pHWE across all variants
p <- create_pretty_hist(dt, aes(x=variant_qc.p_value_hwe), x_label='p(HWE)',
    xlim=c(NA,1), save_figure=FALSE) + scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )
ggsave(paste0(PLOTS, TRANCHE, '_11_pHWE_hist.pdf'), p, width=160, height=90, units='mm')
ggsave(paste0(PLOTS, TRANCHE, '_11_pHWE_hist.jpg'), p, width=160, height=90, units='mm')

# cumulative pHWE
p <- create_pretty_cumulative(dt, aes(x=variant_qc.p_value_hwe), x_label='p(HWE)', threshold=NULL,
    xlim=c(NA,1), save_figure=FALSE) #+ scale_x_log10(
    # breaks = scales::trans_breaks("log10", function(x) 10^x),
    # labels = scales::trans_format("log10", scales::math_format(10^.x))
  # )
ggsave(paste0(PLOTS, TRANCHE, '_11_pHWE_cdf.pdf'), p, width=160, height=90, units='mm')
ggsave(paste0(PLOTS, TRANCHE, '_11_pHWE_cdf.jpg'), p, width=160, height=90, units='mm')

# What other metrics are there?
# Can I also export these other summary metrics?

# For samples

read_and_key <- function(gz_filepath, starts='sample_qc', tmp=FALSE) {
    dt <- fread(cmd = paste('zcat', gz_filepath),
        stringsAsFactors=FALSE, sep='\t', header=TRUE) %>% select(s, starts_with(starts))
    if (tmp) {
        names(dt) <- c('s', paste0('tmp_', names(dt)[-1]))
    }
    setkey(dt, 's')
    return(dt)
}

update_entries <- function(dt) {
    dt <- dt %>% mutate(
        # Counts
        sample_qc.n_called = sample_qc.n_called + tmp_sample_qc.n_called,
        sample_qc.n_not_called = sample_qc.n_not_called + tmp_sample_qc.n_not_called,
        sample_qc.n_filtered = sample_qc.n_filtered + tmp_sample_qc.n_filtered,
        sample_qc.n_hom_ref = sample_qc.n_hom_ref + tmp_sample_qc.n_hom_ref,
        sample_qc.n_het = sample_qc.n_het + tmp_sample_qc.n_het,
        sample_qc.n_hom_var = sample_qc.n_hom_var + tmp_sample_qc.n_hom_var,
        sample_qc.n_non_ref = sample_qc.n_non_ref + tmp_sample_qc.n_non_ref,
        sample_qc.n_singleton = sample_qc.n_singleton + tmp_sample_qc.n_singleton,
        sample_qc.n_snp = sample_qc.n_snp + tmp_sample_qc.n_snp,
        sample_qc.n_insertion = sample_qc.n_insertion + tmp_sample_qc.n_insertion,
        sample_qc.n_deletion = sample_qc.n_deletion + tmp_sample_qc.n_deletion,
        sample_qc.n_transition = sample_qc.n_transition + tmp_sample_qc.n_transition,
        sample_qc.n_transversion = sample_qc.n_transversion + tmp_sample_qc.n_transversion,
        sample_qc.n_star = sample_qc.n_star + tmp_sample_qc.n_star
        ) %>% select(c(s, starts_with('sample_qc')))
    setkey(dt, 's')
    return(dt)
}

add_ratios <- function(dt) {
    dt <- dt %>% mutate(
        # Ratios
        sample_qc.call_rate = sample_qc.n_called / (sample_qc.n_called + sample_qc.n_not_called + sample_qc.n_filtered),
        sample_qc.r_ti_tv = ifelse(sample_qc.n_transversion == 0, NA, sample_qc.n_transition / sample_qc.n_transversion),
        sample_qc.r_het_hom_var = ifelse(sample_qc.n_hom_var == 0, NA, sample_qc.n_het / sample_qc.n_hom_var),
        sample_qc.r_insertion_deletion = ifelse(sample_qc.n_deletion == 0, NA, sample_qc.n_insertion / sample_qc.n_deletion),
    ) %>% select(c('s', starts_with('sample_qc')))
    setkey(dt, 's')
    return(dt)
}

dt <- read_and_key(SAMPLE_QC_FILE)

for (CHR in c(seq(2,22), "X")) {
    # Input files
    cat(paste0("chromosome ", CHR, "\n"))
    SAMPLE_QC_FILE <- paste0(
        "/well/lindgren/UKBIOBANK/dpalmer/wes_", TRANCHE,
        "/ukb_wes_qc/data/samples/11_sample_metrics_for_plotting_chr", CHR, ".tsv.gz")

    dt_tmp <- read_and_key(SAMPLE_QC_FILE, tmp=TRUE)
    dt <- merge(dt, dt_tmp)
    dt <- update_entries(dt)
}

dt <- add_ratios(dt)

# Could include boxplots - I think here I can just include batch.
pdf("test.pdf")
p <- ggplot(dt, aes(y=sample_qc.r_ti_tv)) + geom_boxplot() + theme_classic()
print(p)
dev.off()

# rTiTv
create_pretty_hist(dt, aes(x=sample_qc.r_ti_tv), x_label='Transition/Transversion', save_figure=save_figures,
    file=paste0(PLOTS, TRANCHE, '_11_ti_tv_hist'))

# rHetHomVar
create_pretty_hist(dt, aes(x=sample_qc.r_het_hom_var), x_label='(Het)/(Hom var)', save_figure=save_figures,
    file=paste0(PLOTS, TRANCHE, '_11_het_hom_var_hist'))

# rInsertionDeletion
create_pretty_hist(dt, aes(x=sample_qc.r_insertion_deletion), x_label='Insertion/Deletion', save_figure=save_figures,
    file=paste0(PLOTS, TRANCHE, '_11_insertion_deletion'))

# Number of singletons
create_pretty_hist(dt, aes(x=sample_qc.n_singleton), x_label='Number of singletons', save_figure=save_figures,
    file=paste0(PLOTS, TRANCHE, '_11_n_singleton_hist'))

# Call rate across all samples
create_pretty_hist(dt, aes(x=sample_qc.call_rate), x_label='Call Rate', save_figure=save_figures,
    file=paste0(PLOTS, TRANCHE, '_11_sample_callRate_hist'))

# mean DP all variants
create_pretty_hist(dt, aes(x=sample_qc.dp_stats.mean), x_label='Mean DP', save_figure=save_figures,
    file=paste0(PLOTS, TRANCHE, '_11_sample_mean_DP_hist'))

# mean GQ all variants
create_pretty_hist(dt, aes(x=sample_qc.gq_stats.mean), x_label='Mean GQ', save_figure=save_figures,
    file=paste0(PLOTS, TRANCHE, '_11_sample_mean_GQ_hist'))

# Plot the distribution of these too.
# Compare to the ExAC plots

