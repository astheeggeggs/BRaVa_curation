library(ggplot2)
library(ggsci)
library(dplyr)
library(data.table)

# Plotting functions:
source('utils/pretty_plotting.r')
# Thresholds and plotting file locations defined in r_options.r
source("utils/r_options.r")
source("utils/helpers.r")

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("--tranche", default='200k', help = "Which exome sequencing tranche?")
args <- parser$parse_args()

TRANCHE <- args$tranche

# Inputs
SAMPLE_BEFORE_AFTER_COMBINED_QC_FILE <- paste0('/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE, '/ukb_wes_qc/data/samples/09_final_qc.before.after.combined.samples.tsv'

# Outputs
FINAL_SAMPLE_SUMMARY <- paste0('/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE, '/ukb_wes_qc/data/samples/09_final_sample.summary.tsv')
FINAL_SAMPLE_LIST <- paste0('/well/lindgren/UKBIOBANK/dpalmer/wes_', TRANCHE, '/ukb_wes_qc/data/samples/09_final_qc.keep.BRaVa.sample_list')

dt_after <- fread(SAMPLE_BEFORE_AFTER_COMBINED_QC_FILE, sep='\t', stringsAsFactors=FALSE, header=TRUE) %>%
  filter(phase=='After Variant QC')

dt_keep <- dt_after[, 's']
print(paste0("Started with: ", nrow(dt_keep), " samples"))

# r_ti_tv
dt_keep_ti_tv <- group_by(dt_after, population) %>%
  summarise(
    mean = mean(sample_qc.r_ti_tv),
    sd = sd(sample_qc.r_ti_tv),
    median = median(sample_qc.r_ti_tv),
    mad = mad(sample_qc.r_ti_tv)
  )

# r_het_hom_var
dt_keep_het_hom_var <- group_by(dt_after, population) %>%
  summarise(
    mean = mean(sample_qc.r_het_hom_var),
    sd = sd(sample_qc.r_het_hom_var),
    median = median(sample_qc.r_het_hom_var),
    mad = mad(sample_qc.r_het_hom_var)
  )

# r_insertion_deletion
dt_keep_insertion_deletion <- group_by(dt_after, population) %>%
  summarise(
    mean = mean(sample_qc.r_insertion_deletion),
    sd = sd(sample_qc.r_insertion_deletion),
    median = median(sample_qc.r_insertion_deletion),
    mad = mad(sample_qc.r_insertion_deletion)
  )

# n_singletons
dt_keep_n_singletons <- group_by(dt_after, population) %>%
  summarise(
    mean = mean(sample_qc.n_singleton),
    sd = sd(sample_qc.n_singleton),
    median = median(sample_qc.n_singleton),
    mad = mad(sample_qc.n_singleton)
  )

width <- 160
height <- 90
scaling <- 1

# rTiTv
p <- ggplot(dt_after, aes(x=sample_qc.r_ti_tv, color=factor(population))) + 
  geom_density(stat="density") + theme_classic() + 
  labs(color='1000G label', x='Transition/Transversion') +
  scale_color_d3('category20')
ggsave(paste0(PLOTS, '09_r_ti_tv_density', '.jpg'), p, width=width*scaling, height=height*scaling, units='mm')
ggsave(paste0(PLOTS, '09_r_ti_tv_density', '.pdf'), p, width=width*scaling, height=height*scaling, units='mm')

p <- create_pretty_boxplots(dt_after, aes(y=sample_qc.r_ti_tv, x=factor(population)),
  aes(color=factor(population)), y_label='', x_label='Transition/Transversion')
p <- p + geom_errorbar(data=dt_keep_ti_tv, aes(x=factor(population), ymin=median-n_mads*mad,ymax=median+n_mads*mad,y=median),linetype = 3,width = 0.25)
file <- paste0(PLOTS, '09_r_ti_tv_filter')
ggsave(paste0(file, '.jpg'), p, width=width*scaling, height=height*scaling, units='mm')
ggsave(paste0(file, '.pdf'), p, width=width*scaling, height=height*scaling, units='mm')

# rHetHomVar
p <- ggplot(dt_after, aes(x=sample_qc.r_het_hom_var, color=factor(population))) + 
  geom_density(stat="density") + theme_classic() +
  labs(color='1000G label', x='(Het)/(Hom var)') +
  scale_color_d3('category20')
ggsave(paste0(PLOTS, '09_r_Het_HomVar_density', '.jpg'), p, width=width*scaling, height=height*scaling, units='mm')
ggsave(paste0(PLOTS, '09_r_Het_HomVar_density', '.pdf'), p, width=width*scaling, height=height*scaling, units='mm')

p <- create_pretty_boxplots(dt_after, aes(y=sample_qc.r_het_hom_var, x=factor(population)),
  aes(color=factor(population)), y_label='', x_label='(Het)/(Hom var)')
p <- p + geom_errorbar(data=dt_keep_het_hom_var, aes(x=factor(population), ymin=median-n_mads*mad,ymax=median+n_mads*mad,y=median),linetype = 3,width = 0.25)
file <- paste0(PLOTS, '09_r_Het_HomVar_filter')
ggsave(paste0(file, '.jpg'), p, width=width*scaling, height=height*scaling, units='mm')
ggsave(paste0(file, '.pdf'), p, width=width*scaling, height=height*scaling, units='mm')

# rInsertionDeletion
p <- ggplot(dt_after, aes(x=sample_qc.r_insertion_deletion, color=factor(population))) +
  geom_density(stat="density") + theme_classic() +
  labs(color='1000G label', x='Insertion/Deletion') +
  scale_color_d3('category20')
ggsave(paste0(PLOTS, '09_r_insertion_deletion_density', '.jpg'), p, width=width*scaling, height=height*scaling, units='mm')
ggsave(paste0(PLOTS, '09_r_insertion_deletion_density', '.pdf'), p, width=width*scaling, height=height*scaling, units='mm')

p <- create_pretty_boxplots(dt_after, aes(y=sample_qc.r_insertion_deletion, x=factor(population)),
  aes(color=factor(population)), y_label='', x_label='Insertion/Deletion')
p <- p + geom_errorbar(data=dt_keep_insertion_deletion, aes(x=factor(population), ymin=median-n_mads*mad,ymax=median+n_mads*mad,y=median),linetype = 3,width = 0.25)
file <- paste0(PLOTS, '09_r_insertion_deletion_filter')
ggsave(paste0(file, '.jpg'), p, width=width*scaling, height=height*scaling, units='mm')
ggsave(paste0(file, '.pdf'), p, width=width*scaling, height=height*scaling, units='mm')

# Number of singletons
p <- ggplot(dt_after, aes(x=sample_qc.n_singleton, color=factor(population))) + geom_density(stat="density") + theme_classic() + labs(color='1000G label', x='Number of singletons') 
ggsave(paste0(PLOTS, '09_n_singletons_density', '.jpg'), p, width=width*scaling, height=height*scaling, units='mm')
ggsave(paste0(PLOTS, '09_n_singletons_density', '.pdf'), p, width=width*scaling, height=height*scaling, units='mm')

p <- create_pretty_boxplots(dt_after, aes(y=sample_qc.n_singleton, x=factor(population)),
  aes(color=factor(population)), y_label='', x_label='Number of Singletons')
p <- p + geom_errorbar(data=dt_keep_n_singletons, aes(x=factor(population), ymin=0,ymax=median+20*mad,y=median),linetype = 3,width = 0.25)
file <- paste0(PLOTS, '09_n_Singletons_filter')
ggsave(paste0(file, '.jpg'), p, width=width*scaling, height=height*scaling, units='mm')
ggsave(paste0(file, '.pdf'), p, width=width*scaling, height=height*scaling, units='mm')

# r_ti_tv
dt_keep_ti_tv <- dt_keep_ti_tv %>% inner_join(dt_after, by='population') %>% 
  filter(((sample_qc.r_ti_tv >= (median - n_mads*mad)) & 
     (sample_qc.r_ti_tv <= (median + n_mads*mad))) | is.na(mad))

# r_het_hom_var
dt_keep_het_hom_var <- dt_keep_het_hom_var %>% inner_join(dt_after, by='population') %>% 
  filter(((sample_qc.r_het_hom_var >= (median - n_mads*mad)) & 
     (sample_qc.r_het_hom_var <= (median + n_mads*mad))) | is.na(mad))

# r_insertion_deletion
dt_keep_insertion_deletion <- dt_keep_insertion_deletion %>% inner_join(dt_after, by='population') %>% 
  filter(((sample_qc.r_insertion_deletion >= (median - n_mads*mad)) & 
     (sample_qc.r_insertion_deletion <= (median + n_mads*mad))) | is.na(mad))

# n_singletons
dt_keep_n_singletons <- dt_keep_n_singletons %>% inner_join(dt_after, by='population') %>% 
  filter((sample_qc.n_singleton <= (median + n_mads_singleton*mad)) | is.na(mad))

dt_keep <- dt_keep_ti_tv %>% inner_join(dt_after, by='s') %>% select(s)
print(paste0("Remove Ti/Tv outliers: ", nrow(dt_keep), " samples remain"))

dt_keep <- dt_keep_het_hom_var %>% inner_join(dt_keep, by='s') %>% select(s)
print(paste0("Remove Het/HomVar outliers: ", nrow(dt_keep), " samples remain"))

dt_keep <- dt_keep_insertion_deletion %>% inner_join(dt_keep, by='s') %>% select(s)
print(paste0("Remove Ins/Del outliers: ", nrow(dt_keep), " samples remain"))

dt_keep <- dt_keep_n_singletons %>% inner_join(dt_keep, by='s') %>% select(s)
print(paste0("Remove n_singletons outliers: ", nrow(dt_keep), " samples remain"))

summary_fun <- function(dt) {
  return(
    dt %>% group_by(population) %>%
    summarise("Samples"=n())
    )
}

dt_after_summary <- summary_fun(dt_after)
dt_final_sample_summary <- rbind(
  summary_fun(dt_after)$Samples,
  dt_after_summary$Samples - summary_fun(dt_keep_ti_tv)$Samples,
  dt_after_summary$Samples - summary_fun(dt_keep_het_hom_var)$Samples,
  dt_after_summary$Samples - summary_fun(dt_keep_insertion_deletion)$Samples,
  dt_after_summary$Samples - summary_fun(dt_keep_n_singletons)$Samples,
  summary_fun(merge(dt_keep, dt_after))$Samples
  )

dt_final_sample_summary <- cbind(
  data.table(
    Filter = c(
      "Samples after population filters",
      paste0("Within batch Ti/Tv ratio outside ", n_mads, " standard deviations"),
      paste0("Within batch Het/HomVar ratio outside ", n_mads, " standard deviations"),
      paste0("Within batch Insertion/Deletion ratio outside ", n_mads, " standard deviations"),
      paste0("n singletons > ", n_mads_singleton , " median absolute deviations"),
      "Samples after final sample filters"
      )
    ),
  dt_final_sample_summary)

names(dt_final_sample_summary) <- c("Filter", summary_fun(dt_after)$population)

fwrite(dt_final_sample_summary, file=FINAL_SAMPLE_SUMMARY, quote=FALSE, row.names=FALSE, col.names=TRUE, sep='\t')

# write out
fwrite(dt_keep, file=FINAL_SAMPLE_LIST, quote=FALSE, row.names=FALSE, col.names=FALSE)

