library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(latex2exp)

source("~/Repositories/BRaVa_curation/meta_analysis_utils.r")
system("mkdir -p ~/Repositories/BRaVa_curation/data/meta_analysis")

data_directory <- "~/Repositories/BRaVa_curation/data"
setwd(data_directory)

chr <- c(seq(1,22), "X")
pops <- c("AFR", "AMR", "EAS", "EUR", "SAS")
case_control_threshold <- 100
all_phenotypes <- FALSE
tests <- c("Burden", "SKAT", "SKAT-O")

if (all_phenotypes) {
    # Determine the available phenotypes
    outputs <- system("dx ls brava/outputs/step2", intern=TRUE)
    phenotypes <- unique(gsub(paste0("^chr[0-9X]+_(.*)_[A-Z]{3}_?F?M?.txt.*"), "\\1", outputs))
    if (any(grepl("/", phenotypes))) {
        phenotypes <- phenotypes[-grep("/", phenotypes)]
    }
} else {
    phenotypes <- "Coronary_artery_disease"
}

for (p in phenotypes)
{
    cat(paste0(p, "...\nchromosome "))
    dt_list <- list()
    dt_list_cauchy <- list()
    for (c in chr)
    {
        files <- paste0("chr", c, "_", p, "_", pops)
        cat(paste0(c, "..."))
        for (file in files) {
            file_info <- search_for_files(file)
            if (is.null(file_info)) {
                next
            }
            dt_list[[file]]<- fread(cmd = ifelse(file_info$gz,
                paste0("gzcat ", file_info$gene_file),
                paste0("cat ", file_info$gene_file)))
            dt_list[[file]]$study <- file
            dt_tmp <- add_N(file_info, dt_list[[file]])
            dt_list_cauchy[[file]] <- dt_tmp$dt_gene %>% filter(Group == "Cauchy")
            dt_list[[file]] <- dt_tmp$dt_gene %>% filter(Group != "Cauchy")
        }
    }

    binary <- dt_tmp$binary
    dt <- rbindlist(dt_list)
    dt <- dt %>% filter((!is.na(Pvalue)) & (!is.na(Pvalue_SKAT)) & (!is.na(Pvalue_Burden)))
    
    if (binary) {
        dt <- dt %>% filter(
            N_case > case_control_threshold,
            N_control > case_control_threshold)
    }

    dt_meta <- list()
    for (test in tests)
    {
        Pvalue_col <- ifelse(test == "SKAT-O", "Pvalue", paste0("Pvalue_", test))        
        dt_meta[[test]] <- run_cauchy(
            dt %>% group_by(Region, study),
            "N_eff", "Stat", Pvalue_col, "Pvalue") %>% 
        mutate(type="Cauchy")
    }

    dt_check <- rbindlist(dt_list_cauchy)
    setkeyv(dt_check, c("Region", "study"))
    for (test in tests) {
        dt_meta[[test]] <- data.table(dt_meta[[test]])
        setkeyv(dt_meta[[test]], c("Region", "study"))
    }

    dt_plot <- merge(dt_meta[["Burden"]], dt_check)
    dt_plot <- dt_plot %>% mutate(study = gsub(".*_", "", study))
    pl <- ggplot(dt_plot, aes(x=-log10(Pvalue.x), y=-log10(Pvalue_Burden))) +

    geom_bin_2d(bins=100) + geom_abline(intercept=0,slope=1,col='red',lwd=0.3) +
    theme_bw() +
    labs(
        x=TeX("$-\\log_{10}(P_{Cauchy})$"),
        y=TeX("$-\\log_{10}(P_{SAIGE\\; Cauchy})$")
        ) + facet_wrap(~study)
    print(pl)
}
