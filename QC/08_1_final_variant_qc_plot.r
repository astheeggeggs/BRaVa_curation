library(ggplot2)
library(ggsci)
library(data.table)

source("utils/r_options.r")
source("utils/pretty_plotting.r")

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("--variant_qc_file", required=TRUE, help="Path to VARIANT_QC_FILE prefix to loop over chromosomes")
args <- parser$parse_args()

VARIANT_QC_FILE <- args$variant_qc_file

dt_list <- list()

for (pop in c("AFR", "AMR", "EAS", "EUR", "SAS")) {
    dt_list[[pop]] <- fread(cmd = paste('zcat', COMBINED_VARIANT_QC_FILE), header=TRUE, sep='\t')
    dt_list[[pop]] <- fread(COMBINED_VARIANT_QC_FILE)
    dt_list[[pop]]$pop <- pop
}

dt <- rbindlist(dt_list)
fwrite(dt, file=paste0(VARIANT_QC_FILE, "_combined.tsv"), sep='\t', quote=FALSE)
    
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

