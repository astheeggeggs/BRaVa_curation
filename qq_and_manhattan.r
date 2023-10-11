library(data.table)
library(dplyr)
library(latex2exp)
library(stringr)

source("QC/utils/pretty_plotting.r")

read_and_create_qq_gene <- function(
	phenotype, outdir, indir, input_regexp, population, sex,
	by_chr=FALSE, save=TRUE, gene_mapping="data/gene_mapping.txt.gz",
	protein_coding="data/protein_coding_genes.txt.gz", filter_to_protein_coding=TRUE)
{
	# Mapping file to obtain gene-start positions for each gene.
	# The default mapping file is taken from Biomart, build 38 with stable ID, start and end, chr, and gene name.
	if (grepl("*.gz$", gene_mapping)) {
		dt_gene <- fread(cmd = paste("gzcat", gene_mapping))
	} else {
		dt_gene <- fread(gene_mapping)
	}
	names(dt_gene) <- c("GeneID", "GeneIDversion", "start", "stop", "chr", "GeneName_biomart")

	if (filter_to_protein_coding) {
		# Filter to protein coding
		if (grepl("*.gz$", protein_coding)) {
			dt_gene_protein <- fread(cmd = paste("gzcat", protein_coding))
		} else {
			dt_gene_protein <- fread(protein_coding)
		}
		dt_gene_protein <- dt_gene_protein %>% filter(chromosome_name %in% c(seq(1,22), "X")) %>% select(-hgnc_id)
		dt_gene_protein[, GeneID:=ensembl_gene_id]
		dt_gene_protein[, ensembl_gene_id:=NULL]
		setkey(dt_gene_protein, "GeneID")
		dt_gene_protein <- dt_gene_protein %>% filter(gene_biotype=="protein_coding")
		setkey(dt_gene, "GeneID")
		dt_gene <- merge(dt_gene, dt_gene_protein)
	}
	
	setkeyv(dt_gene, c("GeneID", "chr"))

	outdir <- ifelse(grepl("\\/$", outdir), gsub("\\/.*$", "\\/", outdir), paste0(outdir, "/"))
	files <- dir(indir, pattern=input_regexp, full.names=TRUE)
	print(files)
	dt <- list()
	
	for (file in files) {
		chr <- gsub(".*chr([0-9X]+).*", "\\1", file)
		cat(paste0("chromosome ", gsub(".*chr([0-9X]+).*", "\\1", file), "..."))
		if (grepl("*.gz$", file)) {
			dt[[chr]] <- fread(cmd=paste("gzcat", file))
		} else {
			dt[[chr]] <- fread(file)
		}
		dt[[chr]]$chr <- chr
	}
	
	dt <- rbindlist(dt)
	print(names(dt))
	dt[, Burden := ifelse(is.na(Pvalue_Burden), -log10(Pvalue), -log10(Pvalue_Burden))]
	dt[, SKAT := ifelse(is.na(Pvalue_SKAT), -log10(Pvalue), -log10(Pvalue_SKAT))]
	dt[, `SKAT-O` := -log10(Pvalue)]
	dt[, GeneID := Region]
	dt[, Region := NULL]
	setkeyv(dt, c("GeneID", "chr"))
	dt <- merge(dt, dt_gene)

	ribbon_p <- 0.95

	dt <- melt(
		dt,
		id.vars=c("chr", "GeneID", "Group", "max_MAF", "GeneName_biomart"),
		measure.vars=c("Burden", "SKAT", "SKAT-O"),
		variable.name = "test" , value.name="Pvalue"
		)

	if (by_chr) {
		dt <- data.table(
			dt %>% arrange(chr, test, Group, max_MAF, desc(Pvalue)) %>% select(-c(chr, test, Group, max_MAF)),
			dt %>% group_by(chr, test, Group, max_MAF) %>% arrange(desc(Pvalue)) %>% 
				summarize(
					Pvalue_expected = -log10(seq(1, n())/(n() + 1)),
					clower = -log10(qbeta(p = (1 - ribbon_p) / 2, shape2 = n():1, shape1 = 1:n())),
					cupper = -log10(qbeta(p = (1 + ribbon_p) / 2, shape2 = n():1, shape1 = 1:n()))
				)
			)
	} else {
		dt <- data.table(
			dt %>% arrange(test, Group, max_MAF, desc(Pvalue)) %>% select(-c(test, Group, max_MAF)),
			dt %>% group_by(test, Group, max_MAF) %>% arrange(desc(Pvalue)) %>% 
				summarize(
					Pvalue_expected = -log10(seq(1, n())/(n() + 1)),
					clower = -log10(qbeta(p = (1 - ribbon_p) / 2, shape2 = n():1, shape1 = 1:n())),
					cupper = -log10(qbeta(p = (1 + ribbon_p) / 2, shape2 = n():1, shape1 = 1:n()))
				)
			)
	}

	# First, run everything without Cauchy
	dt_to_plot <- dt %>% filter(Group!="Cauchy")
	max_MAFs <- unique(dt_to_plot$max_MAF)
	variant_classes <- unique(dt_to_plot$Group)
	phenotype_plot <- gsub("_", " ", paste0(toupper(substring(phenotype, 1,1)), substring(phenotype, 2)))
	cat(paste0("\nPhenotype: ", phenotype_plot, "\n"))

	if (save) {
		sex <- ifelse(sex=="both_sexes", "", ifelse(sex=="XX", "_F", "_M"))
		pdf(
			file=paste0(outdir, phenotype, "_gene_",
				ifelse(by_chr, "chr_", ""), population, sex, "_qq.pdf"),
			width=ifelse(by_chr, 15, 4.5),
			height=ifelse(by_chr, 10, 3.5)
			)
	}

	for (c in variant_classes) {
		for (m in max_MAFs) {
			max_MAF_plot <- ifelse(
				grepl("e", as.character(m)),
				gsub("([0-9\\.]+)e(-)*0*([1-9][0-9]*)", "$\\1\\\\times 10^{\\2\\3}$", as.character(m)),
				paste0("$", as.character(m), "$"))
			variant_class_plot <- gsub("_", " ", gsub("[\\|;]", ", ", c))
			p <- create_pretty_qq_plot(
				plot_title=phenotype_plot,
				plot_subtitle=TeX(paste0(variant_class_plot, "; max MAF = ", max_MAF_plot)),
				cex_labels=2,
				dt_to_plot %>% filter(Group==c, max_MAF==m) %>% mutate(labels=GeneName_biomart), aes(x=Pvalue_expected, y=Pvalue, color=test),
				save_figure=FALSE, n_to_include=10,
				x_label=TeX(r'($-\log_{10}(\mathit{p}_{expected})$)'), 
				y_label=TeX("$-\\log_{10}(\\mathit{p}_{observed})$"),
				key_cols=c("test", "Pvalue"),
				aes_ribbon = aes(ymin=clower, ymax=cupper),
				width=170,
				height=120,
				by_chr=by_chr
			)
		}
	}

	# Finally, create the Cauchy QQ plots
	dt_to_plot <- dt %>% filter(Group=="Cauchy")
	create_pretty_qq_plot(
				plot_title=phenotype_plot,
				plot_subtitle="Cauchy combination",
				cex_labels=2,
				dt_to_plot %>% mutate(labels=GeneName_biomart), aes(x=Pvalue_expected, y=Pvalue, color=test),
				save_figure=FALSE, n_to_include=10,
				x_label=TeX("$-\\log_{10}(\\mathit{p}_{expected})$"), 
				y_label=TeX("$-\\log_{10}(\\mathit{p}_{observed})$"),
				key_cols=c("test", "Pvalue"),
				aes_ribbon = aes(ymin=clower, ymax=cupper),
				width=170,
				height=120,
				by_chr=by_chr
			)
	if (save) { dev.off() }
	gc()
}

read_and_create_qq_variant <- function(
	phenotype, outdir, indir, input_regexp, population, sex,
	save=TRUE, split_by_MAF=FALSE, by_chr=FALSE
	)
{
	outdir <- ifelse(grepl("\\/$", outdir), gsub("\\/.*$", "\\/", outdir), paste0(outdir, "/"))

	# Create the QQ plots
	files <- dir(indir, pattern=input_regexp, full.names=TRUE)
	dt <- list()
	
	for (file in files) {
		chr <- gsub(".*chr([0-9X]+).*", "\\1", file)
		cat(paste0("chromosome ", gsub(".*chr([0-9X]+).*", "\\1", file), "..."))
		if (grepl("*.gz$", file)) {
			dt[[chr]] <- fread(cmd=paste("gzcat", file))
		} else {
			dt[[chr]] <- fread(file)
		}
		dt[[chr]]$chr <- chr
	}
	dt <- rbindlist(dt)
	if (nrow(dt) < 1000) { return() }
	dt[, Pvalue := -log10(as.numeric(`p.value`))]

	if("N" %in% names(dt)) {
		dt <- dt %>% filter(((AC_Allele2 >= 20) | (AC_Allele2 >= (2*N - 20))), CHR != "UR")
	} else {
		dt <- dt %>% filter(((AC_Allele2 >= 20) | (AC_Allele2 >= (2*(N_case+N_ctrl) - 20))), CHR != "UR")
	}

	if (nrow(dt) < 1000) { return() }
	ribbon_p <- 0.95
	dt <- data.table(
		dt %>% arrange(desc(Pvalue)),
		dt %>% arrange(desc(Pvalue)) %>% 
			summarize(
				Pvalue_expected = -log10(seq(1, n())/(n() + 1)),
				clower = -log10(qbeta(p = (1 - ribbon_p) / 2, shape2 = n():1, shape1 = 1:n())),
				cupper = -log10(qbeta(p = (1 + ribbon_p) / 2, shape2 = n():1, shape1 = 1:n()))
			)
		)
	phenotype_plot <- gsub("_", " ", paste0(toupper(substring(phenotype, 1,1)), substring(phenotype, 2)))
	cat(paste0("\nPhenotype: ", phenotype_plot, "\n"))

	if (save) {
		sex <- ifelse(sex=="both_sexes", "", ifelse(sex=="XX", "_F", "_M"))
		pdf(
			file=paste0(outdir, phenotype, "_variant_",
				ifelse(by_chr, "chr_", ""), population, sex, "_qq.pdf"),
			width=ifelse(by_chr, 15, 4.5),
			height=ifelse(by_chr, 10, 3.5)
			)
	}

	dt <- dt %>% mutate(MAF = pmin(AF_Allele2, 1-AF_Allele2)) %>% 
		mutate(cols=as.factor(
			ifelse(MAF <1e-4, "$<1\\times 10^{-4}$",
				ifelse(MAF < 1e-3, "$\\[1\\times 10^{-4}, 0.001)$",
					"$\\[0.001, 0.01)$"
					)
				)
			)
		)

	print(dt[which(is.na(dt$Pvalue)),])
	print(dt)
	dt <- data.table(
		dt %>% arrange(cols, desc(Pvalue)) %>% select(-cols),
		dt %>% group_by(cols) %>% arrange(desc(Pvalue)) %>% summarise(Pvalue_expected_MAF = -log10(seq(1, n())/(n() + 1)))
		)

	create_pretty_qq_plot(
		plot_title=phenotype_plot,
		cex_labels=2,
		dt, aes(x=Pvalue_expected, y=Pvalue),
		save_figure=FALSE,
		x_label=TeX("$-\\log_{10}(\\mathit{p}_{expected})$"), 
		y_label=TeX("$-\\log_{10}(\\mathit{p}_{observed})$"),
		key_cols="Pvalue",
		aes_ribbon = aes(ymin=clower, ymax=cupper),
		width=170, height=120,
		by_chr=by_chr
	)

	p <- create_pretty_qq_plot(
		plot_title=phenotype_plot,
		cex_labels=2,
		dt, aes(x=Pvalue_expected_MAF, y=Pvalue, col=as.factor(cols)),
		save_figure=FALSE,
		x_label=TeX("$-\\log_{10}(\\mathit{p}_{expected})$"), 
		y_label=TeX("$-\\log_{10}(\\mathit{p}_{observed})$"),
		key_cols="Pvalue",
		include_qq_ribbon=FALSE,
		width=170, height=120,
		by_chr=by_chr, print_p=FALSE
	)
	p <- p + scale_color_discrete(labels = TeX(levels(dt$cols)))
	print(p)

	if (save) { dev.off() }
	gc()
}

read_and_create_gene_manhattan <- function(
	phenotype, outdir, indir, input_regexp, population, sex,
	save=TRUE,
	gene_mapping="data/gene_mapping.txt.gz",
	threshold=4,
	significance_T=0.05/20000
	)
{
	# Mapping file to obtain gene-start positions for each gene.
	# The default mapping file is taken from Biomart, build 38 with stable ID, start and end, chr, and gene name.
	if (grepl("*.gz$", gene_mapping)) {
		dt_gene <- fread(cmd = paste("gzcat", gene_mapping))
	} else {
		dt_gene <- fread(gene_mapping)
	}
	names(dt_gene) <- c("GeneID", "GeneIDversion", "start", "stop", "chr", "GeneName_biomart")
	setkeyv(dt_gene, c("GeneID", "chr"))

	outdir <- ifelse(grepl("\\/$", outdir), gsub("\\/.*$", "\\/", outdir), paste0(outdir, "/"))
	files <- dir(indir, pattern=input_regexp, full.names=TRUE)
	dt <- list()
	
	for (file in files) {
		chr <- gsub(".*chr([0-9X]+).*", "\\1", file)
		cat(paste0("chromosome ", gsub(".*chr([0-9X]+).*", "\\1", file), "..."))
		if (grepl("*.gz$", file)) {
			dt[[chr]] <- fread(cmd=paste("gzcat", file))
		} else {
			dt[[chr]] <- fread(file)
		}
		dt[[chr]]$chr <- chr
	}
	dt <- rbindlist(dt)

	dt[, Burden := ifelse(is.na(Pvalue_Burden), -log10(Pvalue), -log10(Pvalue_Burden))]
	dt[, SKAT := ifelse(is.na(Pvalue_SKAT), -log10(Pvalue), -log10(Pvalue_SKAT))]
	dt[, `SKAT-O` := -log10(Pvalue)]
	dt[, GeneID := Region]
	dt[, Region := NULL]
	setkeyv(dt, c("GeneID", "chr"))
	dt <- merge(dt, dt_gene)

	dt <- melt(
		dt,
		id.vars=c("chr", "GeneID", "Group", "max_MAF", "GeneName_biomart", "chr", "start", "stop"),
		measure.vars=c("Burden", "SKAT", "SKAT-O"),
		variable.name = "test" , value.name="Pvalue"
		)

	# First, run everything without Cauchy
	dt_to_plot <- dt %>% filter(Group!="Cauchy")
	max_MAFs <- unique(dt_to_plot$max_MAF)
	variant_classes <- unique(dt_to_plot$Group)
	phenotype_plot <- gsub("_", " ", paste0(toupper(substring(phenotype, 1,1)), substring(phenotype, 2)))
	cat(paste0("\nPhenotype: ", phenotype_plot, "\n"))

	if (save) {
		sex <- ifelse(sex=="both_sexes", "", ifelse(sex=="XX", "_F", "_M"))
		pdf(
			file=paste0(outdir, phenotype, "_gene_", population, sex, "_manhattan.pdf"),
			width=12,
			height=3.5
			)
	}

	for (c in variant_classes) {
		for (m in max_MAFs) {
			max_MAF_plot <- ifelse(
				grepl("e", as.character(m)),
				gsub("([0-9\\.]+)e(-)*0*([1-9][0-9]*)", "$\\1\\\\times 10^{\\2\\3}$", as.character(m)),
				paste0("$", as.character(m), "$"))
			variant_class_plot <- gsub("_", " ", gsub("[\\|;]", ", ", c))
			for (t in unique(dt_to_plot$test)) {
				dt_tmp <- dt_to_plot %>% filter(test==t, max_MAF==m, Group==c) %>% select(chr, start, Pvalue, GeneName_biomart)
				make_manhattan_plot(
					dt_tmp$chr, dt_tmp$start, dt_tmp$Pvalue, label=dt_tmp$GeneName_biomart,
					title=phenotype_plot,
					subtitle=TeX(paste0(t, "; ", variant_class_plot, "; max MAF = ", max_MAF_plot)),
					threshold=threshold,
					significance_T=significance_T,
					save_figure=FALSE,
					print_p=TRUE,
					log_p_vals=TRUE
				)
			}
		}
	}

	# Also, need to include the Cauchy Manhattan plot
	dt_to_plot <- dt %>% filter(Group=="Cauchy")
	for (t in unique(dt_to_plot$test)) {
		dt_tmp <- dt_to_plot %>% filter(test==t) %>% select(chr, start, Pvalue, GeneName_biomart)
		make_manhattan_plot(
			dt_tmp$chr, dt_tmp$start, dt_tmp$Pvalue, label=dt_tmp$GeneName_biomart,
			title=phenotype_plot,
			subtitle="Cauchy combination",
			threshold=threshold,
			significance_T=significance_T,
			save_figure=FALSE,
			print_p=TRUE,
			log_p_vals=TRUE
		)
	}

	if (save) { dev.off() }
	gc()
}

read_and_create_variant_manhattan <- function(
	phenotype, outdir, indir, input_regexp, population, sex,
	save=TRUE, significance_T=0.05/1e6
	)
{
	outdir <- ifelse(grepl("\\/$", outdir), gsub("\\/.*$", "\\/", outdir), paste0(outdir, "/"))
	files <- dir(indir, pattern=input_regexp, full.names=TRUE)
	dt <- list()
	
	for (file in files) {
		chr <- gsub(".*chr([0-9X]+).*", "\\1", file)
		cat(paste0("chromosome ", gsub(".*chr([0-9X]+).*", "\\1", file), "..."))
		if (grepl("*.gz$", file)) {
			dt[[chr]] <- fread(cmd=paste("gzcat", file))
		} else {
			dt[[chr]] <- fread(file)
		}
		dt[[chr]]$chr <- chr
	}
	dt <- rbindlist(dt)
	if (nrow(dt) < 1000) { return() }

	if("N" %in% names(dt)) {
		dt <- dt %>% filter(((AC_Allele2 >= 20) | (AC_Allele2 >= (2*N - 20))), CHR != "UR")
	} else {
		dt <- dt %>% filter(((AC_Allele2 >= 20) | (AC_Allele2 >= (2*(N_case+N_ctrl) - 20))), CHR != "UR")
	}
	phenotype_plot <- gsub("_", " ", paste0(toupper(substring(phenotype, 1,1)), substring(phenotype, 2)))
	cat(paste0("\nPhenotype: ", phenotype_plot, "\n"))
	if (nrow(dt) < 1000) { return() }

	if (save) {
		sex <- ifelse(sex=="both_sexes", "", ifelse(sex=="XX", "_F", "_M"))
		pdf(
			file=paste0(outdir, phenotype, "_variant_", population, sex, "_manhattan.pdf"),
			width=12,
			height=3.5
			)
	}

	make_manhattan_plot(
		contigs=dt$CHR, positions=dt$POS, pvals=as.numeric(dt$`p.value`), label=NULL,
		title=phenotype_plot,
		significance_T=significance_T,
		save_figure=FALSE,
		print_p=TRUE
	)
	if (save) { dev.off() }
	gc()
}


create_brava_qq_and_manhattan <- function(
	outdir="plots", indir="data", save=FALSE,
	wait_for_completion=TRUE, gene_mapping="data/gene_mapping.txt.gz",
	include_by_chr=FALSE, overwrite=FALSE,
	RAP_outputs_folder="brava/outputs/step2/sept2023")
{
	# Determine the available phenotypes
	outputs <- system(paste("dx ls", RAP_outputs_folder), intern=TRUE)
	outdir <- ifelse(grepl("\\/$", outdir), gsub("\\/.*$", "\\/", outdir), paste0(outdir, "/"))
	indir <- ifelse(grepl("\\/$", indir), gsub("\\/.*$", "\\/", indir), paste0(indir, "/"))

	dt_results <- data.table(
		output_file=outputs,
		chromosome=gsub("^(chr[0-9X]+)_.*", "\\1", outputs),
		population=gsub(".*_([A-Z]{3}).*.txt.gz", "\\1", outputs),
		sex=ifelse(grepl("F.txt", outputs), "XX",
				ifelse(grepl("M.txt", outputs), "XY", "both_sexes")),
		group_test=ifelse(grepl("singleAssoc", outputs), FALSE, TRUE),
		phenotype=gsub(paste0("^chr[0-9X]+_(.*)_[A-Z]{3}_?F?M?.txt.*"), "\\1", outputs)
	)

	dt_results <- dt_results %>% filter(!grepl("^chr", phenotype))
	system(paste("mkdir -p", paste0(outdir, "gene")))
	system(paste("mkdir -p", paste0(outdir, "variant")))
	system(paste("mkdir -p", paste0(indir, "gene")))
	system(paste("mkdir -p", paste0(indir, "variant")))

	for (s in unique(dt_results$sex)) {
		print(s)
		dt_results_sex <- dt_results %>% filter(sex==s)
		for (p in unique(dt_results_sex$phenotype)) {
			cat(p, "\n")
			# Determine if the phenotype has all of the chromosomes
			for (pop in unique((dt_results_sex %>% filter(phenotype==p))$pop)) {
				cat(paste0(pop, ": "))
				dt_results_sex_tmp <- dt_results_sex %>% filter(phenotype==p, population==pop)
				if (
					all(paste0("chr", c(seq(1,22), "X")) %in% unique(dt_results_sex_tmp$chromosome))
					| !wait_for_completion ) {
					cat(paste0("All chromosomes present for phenotype...\n",
					p, ", population: ", pop, ", sex: ", s, "\n"))

					# Downloading the data if required
					variant_files <- grep("singleAssoc.txt", dt_results_sex_tmp$output_file, value=TRUE)
					gene_files <- setdiff(dt_results_sex_tmp$output_file, variant_files)
					if (!((
							all(file.exists(paste0(indir, "variant/", variant_files))) & 
							all(file.exists(paste0(indir, "gene/", gene_files)))
						) | (
							all(file.exists(paste0(indir, "variant/", variant_files, ".gz"))) & 
							all(file.exists(paste0(indir, "gene/", gene_files, ".gz")))
						))
					)
					{
						cat(paste0("Downloading files for phenotype: ", p, "\n"))
						files_to_download <- dt_results_sex_tmp$output_file
						files_to_download <- unique(gsub("(chr[0-9,X]+)_*", "chr\\*", files_to_download))
						for (files in files_to_download) {
							if (grepl("singleAssoc", files)) {
								system(paste0("dx download ", RAP_outputs_folder, "/", files, " --output ", indir, "variant/"))
							} else {
								system(paste0("dx download ", RAP_outputs_folder, "/", files, " --output ", indir, "gene/"))
							}
						}
					}

					# Remove whitespace in the phenotype naming
					p_trim <- gsub("_$", "", str_trim(gsub("[[:space:]_]+", "\\_", p)))
					# Run gene level plotting
					if (overwrite | !file.exists(
						paste0(outdir, "gene/", p_trim, "_gene_", pop,
							ifelse(s=="both_sexes", "",
								ifelse(s=="XX", "_F", "_M")), "_qq.pdf"))) {
						cat("Creating gene-level QQ plots\n\n")
						read_and_create_qq_gene(
							phenotype=p_trim, outdir=paste0(outdir, "gene"), indir=paste0(indir, "gene"),
							input_regexp=paste0("chr[0-9X]+_", p, "_", pop, ifelse(s=="both_sexes", "", ifelse(s=="XX", "_F", "_M")), ".txt"),
							pop, s, save=save)
						# Run gene level Manhattan plot
						cat("Creating gene-level Manhattan plots\n\n")
						read_and_create_gene_manhattan(
							phenotype=p_trim, outdir=paste0(outdir, "gene"), indir=paste0(indir, "gene"),
							input_regexp=paste0("chr[0-9X]+_", p, "_", pop, ifelse(s=="both_sexes", "", ifelse(s=="XX", "_F", "_M")), ".txt"),
							pop, s, save=save,
							gene_mapping=gene_mapping, threshold=4, significance_T=0.05/20000)
					} else {
						cat("Gene based plots already present!\n")
					}
					# Run variant level plotting
					cat("Creating variant-level QQ plots\n")
					if (overwrite | !file.exists(
					paste0(outdir, "variant/", p_trim, "_variant_", pop,
							ifelse(s=="both_sexes", "",
								ifelse(s=="XX", "_F", "_M")), "_qq.pdf"))) {
					read_and_create_qq_variant(
						phenotype=p_trim, outdir=paste0(outdir, "variant"), indir=paste0(indir, "variant"),
						input_regexp=paste0("chr[0-9X]+_", p, "_", pop, ifelse(s=="both_sexes", "", ifelse(s=="XX", "_F", "_M")), ".txt"),
						pop, s, save=save)
					if (include_by_chr) {
						# Run gene level plotting split by chr
						cat("Creating gene-level QQ plots, split by chromosome\n\n")
						read_and_create_qq_gene(
							phenotype=p_trim, outdir=paste0(outdir, "gene"), indir=paste0(indir, "gene"),
							input_regexp=paste0("chr[0-9X]+_", p, "_", pop, ifelse(s=="both_sexes", "", ifelse(s=="XX", "_F", "_M")), ".txt"),
							pop, s, save=save, by_chr=TRUE)
						read_and_create_qq_variant(
							phenotype=p_trim, outdir=paste0(outdir, "variant"), indir=paste0(indir, "variant"),
							input_regexp=paste0("chr[0-9X]+_", p, "_", pop, ifelse(s=="both_sexes", "", ifelse(s=="XX", "_F", "_M")), ".txt"),
							pop, s, save=save, by_chr=TRUE)
					} 
					# Run variant level Manhattan plot
					cat("Creating variant-level Manhattan plots\n")
					read_and_create_variant_manhattan(
						phenotype=p_trim, outdir=paste0(outdir, "variant"), indir=paste0(indir, "variant"),
						input_regexp=paste0("chr[0-9X]+_", p, "_", pop, ifelse(s=="both_sexes", "", ifelse(s=="XX", "_F", "_M")), ".txt"),
						pop, s, save=save, significance_T=0.05/1e6
						)
					} else {
						cat("Variant based plots already present!\n")
					}
				} else {
					cat(paste0("Not all chromosomes present for phenotype...\n",
					p, ", population: ", pop, "\n"))
					cat(paste(paste(setdiff(paste0("chr", c(seq(1,22), "X")), unique(dt_results_sex_tmp$chromosome)), collapse=", "), "are missing\n")) 
				}
				for (file in paste0(indir, "variant/", variant_files)) {
					if (grepl(".gz$", file)) {
						system(paste("gzip", file))
					}
				}
				for (file in paste0(indir, "gene/", gene_files)) {
					if (grepl(".gz$", file)) {
						system(paste("gzip", file))
					}
				}
			}
		}
	}
}

create_brava_qq_and_manhattan(save=TRUE, wait_for_completion=FALSE)

# Which traits are missing?
# Run traits that are missing
