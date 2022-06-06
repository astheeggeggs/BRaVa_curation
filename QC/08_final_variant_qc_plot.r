library(ggplot2)
library(ggsci)
library(data.table)

source("utils/r_options.r")
source("utils/pretty_plotting.r")

save_figures <- FALSE

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("--variant_qc_file", required=TRUE, help="Path to VARIANT_QC_FILE prefix to loop over chromosomes")
args <- parser$parse_args()

args <- parser$parse_args()

VARIANT_QC_FILE <- args$variant_qc_file
# VARIANT_QC_FILE <- "/well/lindgren/UKBIOBANK/dpalmer/wes_200k/ukb_wes_qc/data/variants/08_final_qc.variants"

for (pop in c("AFR", "AMR", "EAS", "EUR", "SAS"))
{
    CHR <- 1
    save_figures <- TRUE

    # Output files
    COMBINED_VARIANT_QC_FILE <- paste0(VARIANT_QC_FILE, "_", pop, ".tsv")
    VARIANT_QC_FILE_CHR <- paste0(VARIANT_QC_FILE, "_chr", CHR, "_", pop, ".tsv.bgz")
    dt <- fread(cmd = paste('zcat', VARIANT_QC_FILE_CHR), header=TRUE, sep='\t')
    dt_list <- list()
    dt_list[[1]] <- dt

    for (CHR in c(seq(2,22), "X")) {
        # Input files
        cat(paste0("chromosome ", CHR, "\n"))
        VARIANT_QC_FILE_CHR <- paste0(VARIANT_QC_FILE, "_chr", CHR, "_", pop, ".tsv.bgz")
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
    print(nrow(dt_list[[pop]])
    dt_list[[pop]]$pop <- pop
}
dt <- rbindlist(dt_list)
fwrite(dt, file=paste0(VARIANT_QC_FILE, "_combined.tsv"), sep='\t', quote=FALSE))
    
dt <- fread(paste0(VARIANT_QC_FILE, "_combined.tsv"))
p <- create_pretty_cumulative(dt, aes(x=variant_qc.call_rate, col=pop), x_label="Call Rate", threshold=NULL,
    key_label='', xlim=c(0.9,1), save_figure=FALSE)
ggsave(paste0(PLOTS, '08_callRate_cdf', '.jpg'), p, width=160, height=90, units='mm')
ggsave(paste0(PLOTS, '08_callRate_cdf', '.pdf'), p, width=160, height=90, units='mm')

p <- create_pretty_cumulative(
    dt, aes(x=variant_qc.p_value_hwe, col=pop),
    x_label="p-HWE", threshold=NULL,
    key_label='', xlim=c(1e-20,1), save_figure=FALSE) + 
    scale_x_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
            )
ggsave(paste0(PLOTS, '08_pHWE_cdf', '.jpg'), p, width=160, height=90, units='mm')
ggsave(paste0(PLOTS, '08_pHWE_cdf', '.pdf'), p, width=160, height=90, units='mm')

