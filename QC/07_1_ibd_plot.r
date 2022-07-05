library(ggplot2)
library(ggsci)
library(dplyr)
library(latex2exp)
library(data.table)

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("--ibd_file", required=TRUE, help="Path to IBD_OUTPUT from 07_0_ibd.py")
parser$add_argument("--ibd_threshold", default=0.2, help="PI HAT threshold used for relatedness")
args <- parser$parse_args()

IBD_FILE <- args$ibd_file
IBD_THRESHOLD <- args$ibd_threshold

df <- fread(IBD_FILE, sep='\t', stringsAsFactors=TRUE, header=TRUE, data.table=FALSE)

df$inferred_relationship <- 'Siblings'
df$inferred_relationship[df$ibd.PI_HAT <= IBD_THRESHOLD] <- 'Unrelated'
df$inferred_relationship[df$ibd.Z0 < 0.05 & df$ibd.Z1 < 0.05] <- 'Duplicate/Monozygotic twins'
df$inferred_relationship[df$ibd.Z0 < 0.05 & df$ibd.Z1 > 0.9] <- 'Parent-Offspring'

p <- ggplot(df, aes(x=ibd.Z0)) +
  geom_point(aes(y=ibd.Z1, color=inferred_relationship), size=0.5) +
  geom_abline(intercept=(2-2*IBD_THRESHOLD), slope=-2, linetype='dashed') +
  scale_color_d3('category20', limits=c('Unrelated', 'Siblings', 'Parent-Offspring', 'Duplicate/Monozygotic twins')) +
  labs(x='Proportion of loci with 0 shared alleles',
       y='Proportion of loci with 1 shared allele',
       title="",
       color='Inferred Relationship') +
  scale_x_continuous(breaks=scales::pretty_breaks(n=10)) +
  scale_y_continuous(breaks=scales::pretty_breaks(n=10)) +
  theme_minimal() +
  theme(axis.title.x = element_text(margin=margin(t=10)),
        axis.title.y = element_text(margin=margin(r=10)),
        plot.title = element_text(hjust=0.5))
print(p)
ggsave(paste0(PLOTS, '07_IBD_plot', '.jpg'), p, width=160, height=90, units='mm')
ggsave(paste0(PLOTS, '07_IBD_plot', '.pdf'), p, width=160, height=90, units='mm')
