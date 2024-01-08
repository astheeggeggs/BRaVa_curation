#!/usr/bin/env Rscript
library(data.table)
library(dplyr)
library(ggplot2)
library(argparse)

color_amr <- '#ED1E24'
color_eur <- '#6AA5CD'
color_afr <- '#941494'
color_sas <- '#FF9912'
color_eas <- '#108C44'
color_oth <- '#ABB9B9'
color_mde <- '#33CC33'
color_asj <- 'coral'
color_nfe <- color_eur
color_fin <- '#002F6C'

pop_colors <- c('AFR' = color_afr,
               'AMR' = color_amr,
               'EAS' = color_eas,
               'FIN' = color_fin,
               'EUR' = color_nfe,
               'NEF' = color_nfe,
               'OTH' = color_oth,
               'SAS' = color_sas,
               'MDE' = color_mde,
               'ASJ' = color_asj,
               'uniform' = 'pink',
               'consanguineous' = 'pink',
               'SAS_non_consang' = 'orange')

pop_names <- c('OTH' = 'Other',
              'AFR' = 'African/African-American',
              'AMR' = 'Latino',
              'EAS' = 'East Asian',
              'FIN' = 'Finnish',
              'EUR' = 'European',
              'NFE' = 'European',
              'SAS' = 'South Asian',
              'MDE' = 'Middle Eastern',
              'ASJ' = 'Ashkenazi Jewish',
              'uniform' = 'Uniform',
              'SAS_non_consang' = 'South Asian (F < 0.05)',
              'consanguineous' = 'South Asian (F > 0.05)')

annotation_names <- c(
	"damaging_missense_or_protein_altering" = "Damaging missense or protein altering",
	"non_coding" = "Non-coding",
	"other_missense_or_protein_altering" = "Other missense or protein altering",
	"pLoF" = "pLoF",
	"synonymous" = "Synonymous")

annotation_names <- c(
	"damaging_missense_or_protein_altering" = "Damaging missense or protein altering",
	"damaging_missense" = "Damaging missense or protein altering",
	"non_coding" = "Non-coding",
	"other_missense_or_protein_altering" = "Other missense or protein altering",
	"other_missense" = "Other missense or protein altering",
	"pLoF" = "pLoF",
	"synonymous" = "Synonymous")

# Download the required summary data, looping over chromosomes.
# Sum them up and create barplots of the results
# system("mkdir -p counts")
# system("dx download /Duncan/annotation_counts/* -f -o counts/")
# Mirror the ExAC paper
# Linear and log scale

extract_overall_counts <- function(
	file_grep, count_directory="counts",
	annotation_group=annotation, bin_type="MAF")
{
	files_grep <- gsub("@", "[0-9]*X?", file_grep)
	bin <- bin_type
	dt_list <- list()
	files <- grep(files_grep, dir(count_directory, full.names=TRUE), value=TRUE)
	if (length(files) > 0) {
		for (file in grep(files_grep, dir(count_directory, full.names=TRUE), value=TRUE)) {
			dt_list[[file]] <- fread(file)
		}
		dt <- rbindlist(dt_list)
		annotation_group <- enquo(annotation_group)
		return(dt %>% filter(bin_type==!!bin) %>% 
			group_by(!!annotation_group, bin) %>% 
			summarise(
				chromosome_wide_count = sum(average_count),
				chromosome_wide_variants = sum(variant_count))
			)
	} else {
		return(NULL)
	}
}

levels_list <- list(
	MAF=c("<0.1%", "0.1-1%", "1-5%", ">5%"),
	MAC=c("singletons", "(1,5]",  "(5,10]", "(10,100]" , "(100,1,000]", "(1,000,10,000]", ">10,000"),
	spliceAI=c("<0.2", "[0.2,0.5)", "[0.5,0.8)", ">0.8")
)

main <- function(args) {
	biobank <- args$biobank
	count_dir <- args$count_directory
	out <- args$out
	pops <- c("EAS", "AMR", "AFR", "EUR", "SAS")
	if (!is.null(args$super_population)) {
		if (!((args$super_population) %in% pops)) {
			stop("Error: super population specified is not one of EAS, AMR, AFR, EUR, or SAS.")
		}
		pops <- args$super_population
	}
	# BRaVa annotations
	out_variant_AC <- paste0(gsub(".pdf$", "", out), ".variants_AC_bins.pdf")
	pdf(file=out_variant_AC, width=10, height=6)
	for (annotation_level in c("variant", "transcript")) {
		dt_list <- list()
		for (type in c("MAF", "MAC")) {
			for (pop in pops) {
				dt_tmp <- extract_overall_counts(
					paste0(biobank, ".", pop, ".chr@.BRaVa_annotations_", annotation_level, "_summary.tsv.gz"),
					count_directory=count_dir,
					bin_type=type)
				if (!is.null(dt_tmp)) {
					dt_list[[type]][[pop]] <- dt_tmp
					dt_list[[type]][[pop]]$population <- pop
				}
			}
			dt <- rbindlist(dt_list[[type]])
			dt <- dt %>% filter(annotation != "non_coding")
			
			# Axis treated as discrete variable
			if (is.null(levels_list[[type]])) {
				dt$bin <- factor(dt$bin)
			} else {
				dt$bin <- factor(dt$bin, levels=levels_list[[type]])
			}

			p <- ggplot(data=dt, aes(x=bin, y=chromosome_wide_count, fill=population)) + 
				geom_bar(stat="identity", position=position_dodge()) + 
				facet_wrap(vars(annotation), scales="free",
					labeller=labeller(annotation=annotation_names)) +
				scale_fill_manual(values=pop_colors) + 
				labs(x=paste(type, "bin"), y="Count") + 
				guides(x =  guide_axis(angle = 30)) +
				theme_minimal() + theme(legend.title=element_blank())
			print(p)
			p <- p + scale_y_log10()
			print(p)

			p <- ggplot(data=dt, aes(x=bin, y=chromosome_wide_variants, fill=population)) + 
				geom_bar(stat="identity", position=position_dodge()) + 
				facet_wrap(vars(annotation), scales="free",
					labeller=labeller(annotation=annotation_names)) +
				scale_fill_manual(values=pop_colors) + 
				labs(x=paste(type, "bin"), y="Count") + 
				guides(x =  guide_axis(angle = 30)) +
				theme_minimal() + theme(legend.title=element_blank())
			print(p)
			p <- p + scale_y_log10()
			print(p)
		}
	}
	dev.off()

	# SpliceAI binning
	# just the overall counts
	type <- "spliceAI"
	out_spliceAI <- paste0(gsub(".pdf$", "", out), ".spliceAI_bins.pdf")
	pdf(file=out_spliceAI, width=5, height=3)
	for (annotation_level in c("variant", "transcript"))
	{
		dt_list <- list()
		for (pop in pops) {
			dt_tmp <- dt_list[[type]][[pop]] <- extract_overall_counts(
				paste0(biobank, ".", pop, ".chr@.BRaVa_annotations_", annotation_level, "_summary.tsv.gz"),
				count_directory=count_dir,
				bin_type=type)
			if (!is.null(dt_tmp)) {
				dt_list[[type]][[pop]] <- dt_tmp
				dt_list[[type]][[pop]]$population <- pop
			}
		}
		dt <- rbindlist(dt_list[[type]]) %>% 
			group_by(bin, population) %>% 
			summarise(
				chromosome_wide_count = sum(chromosome_wide_count),
				chromosome_wide_variants = sum(chromosome_wide_variants)
			)
		dt <- dt %>% filter(!(is.na(bin) | bin==""))
		
		# Axis treated as discrete variable
		if (is.null(levels_list[[type]])) {
			dt$bin <- factor(dt$bin)
		} else {
			dt$bin <- factor(dt$bin, levels=levels_list[[type]])
		}

		p <- ggplot(data=dt, aes(x=bin, y=chromosome_wide_count, fill=population)) + 
			geom_bar(stat="identity", position=position_dodge()) + 
			scale_fill_manual(values=pop_colors) + 
			labs(x=paste(type,  "bin"), y="Count") +
			guides(x =  guide_axis(angle = 30)) +
			theme_minimal() + theme(legend.title=element_blank())
		print(p)
		p <- p + scale_y_log10()
		print(p)

		p <- ggplot(data=dt, aes(x=bin, y=chromosome_wide_variants, fill=population)) + 
			geom_bar(stat="identity", position=position_dodge()) + 
			scale_fill_manual(values=pop_colors) + 
			labs(x=paste(type,  "bin"), y="Count") +
			guides(x =  guide_axis(angle = 30)) +
			theme_minimal() + theme(legend.title=element_blank())
		print(p)
		p <- p + scale_y_log10()
		print(p)
	}
	dev.off()
}
# Add arguments
parser <- ArgumentParser()
parser$add_argument("--count_directory", default=NULL, required=TRUE,
    help="Directory containing summary files generated from 'merge_annots_and_AC.r' of the
    form {biobank}.{pop}.chr{chr}.BRaVa_annotations_{variant, transcript}_summary.tsv.gz")
parser$add_argument("--biobank", default=NULL, required=TRUE,
    help="biobank name, note this will only work for files saved of the form
    {biobank}.{pop}.chr{chr}.BRaVa_annotations_{variant, transcript}_summary.tsv.gz")
parser$add_argument("--out", default=NULL, required=TRUE,
    help="Output filepath for pdf file of plots")
parser$add_argument("--spliceAI_bins", default=FALSE, action="store_true")
parser$add_argument("--super_population", default=NULL, required=FALSE)
args <- parser$parse_args()

main(args)

