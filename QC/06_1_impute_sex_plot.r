rm(list=ls())
library(ggplot2)
library(ggsci)
library(dplyr)
library(ggExtra)
library(data.table)

# File locations and plotting locations defined in r_options.r
source("utils/r_options.r")
source("utils/pretty_plotting.r")

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("--impute_sex_table", required=TRUE, help="Path to IMPUTESEX_TABLE from 06_0_impute_sex.py")
parser$add_argument("--y_ncalled", required=TRUE, help="Path to Y_NCALLED from 06_0_impute_sex.py")
parser$add_argument("--sexcheck_list", required=TRUE, help="Path to output file containing sex swaps")
args <- parser$parse_args()

IMPUTESEX_FILE <- args$impute_sex_table
Y_NCALLED_FILE <- args$y_ncalled
SEXCHECK_LIST <- args$sexcheck_list

df <- fread(IMPUTESEX_FILE, sep='\t', stringsAsFactors=FALSE, header=TRUE, data.table=FALSE) %>%
  mutate(imputed_sex=as.factor(ifelse(impute_sex.is_female == TRUE, 'Female', 'Male')))

df_y <- fread(Y_NCALLED_FILE, sep='\t', stringsAsFactors=FALSE, header=TRUE, data.table=FALSE)

df <- merge(df, df_y, by='s')

colors <- pal_d3('category20')(20)[c(1,2)]
fills <- pal_d3('category20')(20)[c(11,12)]

create_pretty_cumulative(df, aes(impute_sex.f_stat), 'F-statistic', T_impute_sex,
    xlim=c(-T_impute_sex0,1.1), title='Cumulative Distribution of F-statistic', save_figure=TRUE, file=paste0(PLOTS,'05_F_stat_cdf'))

p <- ggplot(df, aes(x=impute_sex.f_stat, fill=imputed_sex)) +
  geom_histogram(binwidth=0.025, alpha=0.8, color='#7f7f7f') +
  scale_fill_manual(values=fills, limits=c('Male', 'Female')) +
  labs(x='X chromosome F-statistic',
       y='Count',
       title='',
       fill='Imputed Sex') +
  scale_x_continuous(breaks=scales::pretty_breaks(n=10)) +
  scale_y_continuous(label=scales::comma, breaks=scales::pretty_breaks(n=10)) +
  theme_classic() +
  theme(axis.title.x=element_text(margin=ggplot2::margin(t=10)),
        axis.title.y=element_text(margin=ggplot2::margin(r=10)),
        plot.title=element_text(hjust=0.5)) +
  geom_vline(xintercept=T_impute_sex, linetype='dashed')

print(p)
ggsave(paste0(PLOTS, '06_imputesex_histogram', '.pdf'), p, width=160, height=90, units='mm')
ggsave(paste0(PLOTS, '06_imputesex_histogram', '.jpg'), p, width=160, height=90, units='mm', dpi=500)

if(any(df$phenotype.reported_sex == '' | df$phenotype.reported_sex == 'Not Reported' | df$phenotype.reported_sex == 'Unknown')) {
  df$phenotype.reported_sex[df$phenotype.reported_sex == '' | df$phenotype.reported_sex == 'Not Reported' | df$phenotype.reported_sex == 'Unknown'] <- "Unknown"
}

p <- ggplot(df, aes(x=impute_sex.f_stat, y=phenotype.LOCATION, colour=phenotype.reported_sex)) +
  geom_jitter(width=0, height=0.2, size=1, alpha=0.2, stroke=0.05) + 
  theme_minimal() +
  geom_vline(xintercept=T_impute_sex, linetype='dashed') +
  labs(x='X chromosome F-statistic',
       y='Location',
       color='Reported Sex') 

print(p)
ggsave(paste0(PLOTS, '06_imputesex_scatter_box', '.jpg'), p, width=160, height=90, units='mm')
ggsave(paste0(PLOTS, '06_imputesex_scatter_box', '.pdf'), p, width=160, height=90, units='mm')

df_false <- df %>% filter((impute_sex.f_stat > T_impute_sex & phenotype.reported_sex == 'Female') | (impute_sex.f_stat < T_impute_sex & phenotype.reported_sex == 'Male'))
df_false_plot <- df %>% filter(phenotype.reported_sex == 'Unknown' | (impute_sex.f_stat > T_impute_sex & phenotype.reported_sex == 'Female') | (impute_sex.f_stat < T_impute_sex & phenotype.reported_sex == 'Male'))

# Plots of sex estimates.
p <- ggplot(df, aes(x=impute_sex.f_stat, y=impute_sex.n_called, colour=phenotype.reported_sex)) +
geom_point(size=0.5) + 
labs(x='X chromosome F-statistic', y='Number of calls in Y', color='Reported Sex') +
scale_color_d3('category10') +
scale_x_continuous(breaks=scales::pretty_breaks(n=10)) +
scale_y_continuous(label=scales::comma, breaks=scales::pretty_breaks(n=10)) +
geom_point(data=df_false_plot, aes(x=impute_sex.f_stat, y=impute_sex.n_called), size=0.5) + 
theme_minimal()
print(p)

ggsave(paste0(PLOTS, '06_imputesex_scatter', '.jpg'), p, width=160, height=90, units='mm')
ggsave(paste0(PLOTS, '06_imputesex_scatter', '.pdf'), p, width=160, height=90, units='mm')

df_out <- df_false %>% select(s)
write.table(df_out, file=SEXCHECK_LIST, quote=FALSE, row.names=FALSE, col.names=FALSE)
