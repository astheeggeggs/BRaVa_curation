rm(list=ls())
library(ggplot2)
library(ggsci)
library(dplyr)
library(data.table)

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("--sample_before_prefix", help = "Prefix to the set of files before sample QC, output by 09_final_sample_QC.py")
parser$add_argument("--sample_after_prefix", help = "Prefix to the set of files after sample QC, output by 09_final_sample_QC.py")
parser$add_argument("--output", help = "Output file combining sample metrics before and after curation across chromosomes and populations for plotting")
args <- parser$parse_args()

SAMPLE_BEFORE_QC_FILE <- args$sample_before_prefix
SAMPLE_AFTER_QC_FILE <- args$sample_after_prefix
SAMPLE_BEFORE_AFTER_COMBINED_QC_FILE <- args$output

# Plotting functions:
source('utils/pretty_plotting.r')
# Thresholds and plotting file locations defined in r_options.r
source("utils/r_options.r")
source("utils/helpers.r")

save_figures <- TRUE

CHR <- 1

# Input files
SAMPLE_BEFORE_QC_FILE <-  '/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc_chr@.before.'
SAMPLE_AFTER_QC_FILE <-  '/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc_chr@.after.'

# Output files
SAMPLE_BEFORE_AFTER_COMBINED_QC_FILE <- '/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/samples/09_final_qc.before.after.combined.samples.tsv'

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
    cat(paste0(pop, "...\n"))
    dt_before <- read_and_key(paste0(SAMPLE_BEFORE_QC_FILE , ".", pop, ".samples.tsv.bgz"))
    dt_after <- read_and_key(paste0(SAMPLE_AFTER_QC_FILE), ".", pop, ".samples.tsv.bgz"))

    dt_before <- dt_before %>% mutate(phase = 'Before Variant QC')
    dt_after <- dt_after %>% mutate(phase = 'After Variant QC')

    dt <- bind_rows(dt_before, dt_after) %>%
        mutate(
            phase=factor(phase, levels=c('Before Variant QC', 'After Variant QC')),
            population = pop
            )
    dt$population <- pop
    dt_list[[pop]] <- dt
}
dt <- rbindlist(dt_list)

fwrite(dt, file=SAMPLE_BEFORE_AFTER_COMBINED_QC_FILE, sep='\t')

# Split by 1000G populations.
if (save_figures)
{
    y_labels <- ''
    alpha <- 0.8

    # Number of singletons
    create_pretty_boxplots(dt, aes(y=sample_qc.n_singleton, x=factor(population)),
        aes(color=factor(population)), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='Number of Singletons',
        save_figure=save_figures, file=paste0(PLOTS, '09_nSingletons'), y_label=y_labels,
        alpha=alpha, height=140)

    # rHetHomVar
    create_pretty_boxplots(dt, aes(y=sample_qc.r_het_hom_var, x=factor(population)),
        aes(color=factor(population)), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rHetHomVar',
        save_figure=save_figures, file=paste0(PLOTS, '09_rHetHomVar'), y_label=y_labels,
        alpha=alpha, height=140)

    # rInsertionDeletion
    create_pretty_boxplots(dt, aes(y=sample_qc.r_insertion_deletion, x=factor(population)),
        aes(color=factor(population)), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rInsertionDeletion',
        save_figure=save_figures, file=paste0(PLOTS, '09_rInsertionDeletion'), y_label=y_labels,
        alpha=alpha, height=140)

    # rTiTv
    create_pretty_boxplots(dt, aes(y=sample_qc.r_ti_tv, x=factor(population)),
        aes(color=factor(population)), facet=TRUE, facet_grid=facet_grid(phase~.), x_label='rTiTv',
        save_figure=save_figures, file=paste0(PLOTS, '09_rTiTv_by_centre'), n_ticks=5, y_label=y_labels,
        alpha=alpha, height=140)
}
