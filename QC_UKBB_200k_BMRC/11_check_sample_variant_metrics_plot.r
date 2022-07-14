library(ggplot2)
library(ggsci)
library(data.table)

# Plotting functions:
source('utils/pretty_plotting.r')
# Thresholds and plotting file locations defined in r_options.r
source("utils/r_options.r")

suppressPackageStartupMessages(library("argparse"))

TRANCHE <- '200k'
save_figures <- TRUE

# Input files
SAMPLE_QC_FILE <- paste0("/well/lindgren/UKBIOBANK/dpalmer/wes_", TRANCHE, "/ukb_wes_qc/data/11_sample_metrics_for_plotting")
SAMPLE_QC_TARGET_FILE <- paste0("/well/lindgren/UKBIOBANK/dpalmer/wes_", TRANCHE, "/ukb_wes_qc/data/11_sample_target_interval_metrics_for_plotting")
VARIANT_QC_FILE <- paste0("/well/lindgren/UKBIOBANK/dpalmer/wes_", TRANCHE, "/ukb_wes_qc/data/11_variant_metrics_for_plotting")

# Output files
COMBINED_SAMPLE_QC_FILE <- paste0('/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE, '/ukb_wes_qc/data/11_sample_metrics_for_plotting.tsv')

dt_list <- list()
dt_list[[1]] <- dt

for (pop in c("AFR", "AMR", "EAS", "EUR", "SAS"))
{
    CHR <- 19#1
    # Output files
    COMBINED_VARIANT_QC_FILE <- paste0(VARIANT_QC_FILE, "_", pop, ".tsv")
    VARIANT_QC_FILE_CHR <- paste0(VARIANT_QC_FILE, "_chr", CHR, ".", pop, ".tsv.bgz")
    dt <- fread(cmd = paste('zcat', VARIANT_QC_FILE_CHR), header=TRUE, sep='\t')
    dt_list <- list()
    dt_list[[1]] <- dt

    for (CHR in c(seq(20,20))) {#c(seq(2,22), "X")) {
        # Input files
        cat(paste0("chromosome ", CHR, "\n"))
        VARIANT_QC_FILE_CHR <- paste0(VARIANT_QC_FILE, "_chr", CHR, ".", pop, ".tsv.bgz")
        dt_tmp <- fread(cmd = paste('zcat', VARIANT_QC_FILE_CHR), header=TRUE, sep='\t')
        setkeyv(dt_tmp, c('locus', 'alleles'))
        dt_list[[CHR]] <- dt_tmp
    }

    dt <- rbindlist(dt_list)
    fwrite(dt, file=COMBINED_VARIANT_QC_FILE , sep='\t', quote=FALSE)
}

dt_list <- list()
for (pop in c("AFR", "AMR", "EAS", "EUR", "SAS")) {
    COMBINED_VARIANT_QC_FILE <- paste0(VARIANT_QC_FILE, "_", pop, ".tsv")
    dt_list[[pop]] <- fread(COMBINED_VARIANT_QC_FILE)
    dt_list[[pop]]$pop <- pop
}

dt <- rbindlist(dt_list)
fwrite(dt, file=paste0(VARIANT_QC_FILE, "_combined.tsv"), sep='\t', quote=FALSE)

dt <- fread(paste0(VARIANT_QC_FILE, "_combined.tsv"))

# Variants

# cumulative average call rate
create_pretty_cumulative(dt, aes(x=variant_qc.call_rate, col=pop), x_label="Call Rate", threshold=NULL,
    key_label='', xlim=c(NA,1), save_figure=save_figures, file=paste0(PLOTS, '11_callRate_cdf'))

# pHWE
p <- create_pretty_cumulative(
    dt, aes(x=variant_qc.p_value_hwe, col=pop),
    x_label="p-HWE", threshold=NULL,
    key_label='', xlim=c(1e-20,1), save_figure=FALSE) + 
    scale_x_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
            )
ggsave(paste0(PLOTS, '11_pHWE_cdf', '.jpg'), p, width=160, height=90, units='mm')
ggsave(paste0(PLOTS, '11_pHWE_cdf', '.pdf'), p, width=160, height=90, units='mm')

# cumulative mean GQ
create_pretty_cumulative(dt, aes(x=variant_qc.gq_stats.mean, col=pop), x_label="Mean GQ", threshold=NULL,
    key_label='', save_figure=save_figures, file=paste0(PLOTS, '11_mean_GQ_cdf'))

# cumulative mean DP
create_pretty_cumulative(dt, aes(x=variant_qc.dp_stats.mean, col=pop), x_label="Mean DP", threshold=NULL,
    key_label='', save_figure=save_figures, file=paste0(PLOTS, '11_mean_DP_cdf'))

# Samples

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


dt_list <- list()
for (pop in c("AFR", "AMR", "EAS", "EUR", "SAS"))
{
    CHR <- 19#1
    cat(paste0(pop, "...\n"))
    cat("chromosome 1\n")
    dt <- read_and_key(paste0(SAMPLE_QC_FILE, "_chr", CHR, ".", pop, ".tsv.bgz"))
    dt_target <- read_and_key(paste0(SAMPLE_QC_TARGET_FILE, "_chr", CHR, ".", pop, ".tsv.bgz"))

    for (CHR in seq(20,20)){#c(seq(2,22), "X")) {
        # Input files
        cat(paste0("chromosome ", CHR, "\n"))
        dt_tmp <- read_and_key(paste0(SAMPLE_QC_FILE, "_chr", CHR, ".", pop, ".tsv.bgz"), tmp=TRUE)
        dt_target_tmp <- read_and_key(paste0(SAMPLE_QC_TARGET_FILE, "_chr", CHR, ".", pop, ".tsv.bgz"), tmp=TRUE)

        dt  <- merge(dt, dt_tmp)
        dt_target <- merge(dt_target, dt_target_tmp)

        dt<- update_entries(dt)
        dt_target <- update_entries(dt_target)
    }

    dt <- add_ratios(dt) %>% mutate(phase = 'All QCed Variants')
    dt_target <- add_ratios(dt_target) %>% mutate(phase = 'QCed Variants in Target')
    dt <- bind_rows(dt, dt_target) %>%
        mutate(
            phase=factor(phase, levels=c('All QCed Variants', 'QCed Variants in Target')),
            population = pop
            )
    dt$population <- pop
    dt_list[[pop]] <- dt
}

dt <- rbindlist(dt_list)

fwrite(dt, file=COMBINED_SAMPLE_QC_FILE, sep='\t')
dt <- fread(COMBINED_SAMPLE_QC_FILE)

# Split by 1000G populations.
if (save_figures)
{
    y_labels <- ''
    alpha <- 0.8

    # Number of singletons
    create_pretty_boxplots(dt, aes(y=sample_qc.n_singleton, x=factor(population)),
        aes(color=factor(population)), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='Number of Singletons',
        save_figure=save_figures, file=paste0(PLOTS, '11_nSingletons_by_1000G_label'), y_label=y_labels,
        alpha=alpha, height=140)

    # rHetHomVar
    create_pretty_boxplots(dt, aes(y=sample_qc.r_het_hom_var, x=factor(population)),
        aes(color=factor(population)), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rHetHomVar',
        save_figure=save_figures, file=paste0(PLOTS, '11_rHetHomVar_by_1000G_label'), y_label=y_labels,
        alpha=alpha, height=140)

    # rInsertionDeletion
    create_pretty_boxplots(dt, aes(y=sample_qc.r_insertion_deletion, x=factor(population)),
        aes(color=factor(population)), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rInsertionDeletion',
        save_figure=save_figures, file=paste0(PLOTS, '11_rInsertionDeletion_by_1000G_label'), y_label=y_labels,
        alpha=alpha, height=140)

    # rTiTv
    create_pretty_boxplots(dt, aes(y=sample_qc.r_ti_tv, x=factor(population)),
        aes(color=factor(population)), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rTiTv',
        save_figure=save_figures, file=paste0(PLOTS, '11_rTiTv_by_centre_by_1000G_label'), n_ticks=5, y_label=y_labels,
        alpha=alpha, height=140)
}




# # Could include boxplots - I think here I can just include batch.
# pdf("test.pdf")
# p <- ggplot(dt, aes(y=sample_qc.r_ti_tv)) + geom_boxplot() + theme_classic()
# print(p)
# dev.off()

# # rTiTv
# create_pretty_hist(dt, aes(x=sample_qc.r_ti_tv), x_label='Transition/Transversion', save_figure=save_figures,
#     file=paste0(PLOTS, TRANCHE, '_11_ti_tv_hist'))

# # rHetHomVar
# create_pretty_hist(dt, aes(x=sample_qc.r_het_hom_var), x_label='(Het)/(Hom var)', save_figure=save_figures,
#     file=paste0(PLOTS, TRANCHE, '_11_het_hom_var_hist'))

# # rInsertionDeletion
# create_pretty_hist(dt, aes(x=sample_qc.r_insertion_deletion), x_label='Insertion/Deletion', save_figure=save_figures,
#     file=paste0(PLOTS, TRANCHE, '_11_insertion_deletion'))

# # Number of singletons
# create_pretty_hist(dt, aes(x=sample_qc.n_singleton), x_label='Number of singletons', save_figure=save_figures,
#     file=paste0(PLOTS, TRANCHE, '_11_n_singleton_hist'))

# # Call rate across all samples
# create_pretty_hist(dt, aes(x=sample_qc.call_rate), x_label='Call Rate', save_figure=save_figures,
#     file=paste0(PLOTS, TRANCHE, '_11_sample_callRate_hist'))

# # mean DP all variants
# create_pretty_hist(dt, aes(x=sample_qc.dp_stats.mean), x_label='Mean DP', save_figure=save_figures,
#     file=paste0(PLOTS, TRANCHE, '_11_sample_mean_DP_hist'))

# # mean GQ all variants
# create_pretty_hist(dt, aes(x=sample_qc.gq_stats.mean), x_label='Mean GQ', save_figure=save_figures,
#     file=paste0(PLOTS, TRANCHE, '_11_sample_mean_GQ_hist'))

# # Plot the distribution of these too.
# # Compare to the ExAC plots

