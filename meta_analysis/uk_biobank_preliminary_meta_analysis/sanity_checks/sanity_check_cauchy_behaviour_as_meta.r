library(data.table)
library(dplyr)
library(stringr)
source("~/Repositories/BRaVa_curation/meta_analysis/meta_analysis_utils.r")
system("mkdir -p ~/Repositories/BRaVa_curation/data/meta_analysis")

data_directory <- "~/Repositories/BRaVa_curation/data"
setwd(data_directory)

chr <- c(seq(1,22), "X")
pops <- c("EUR", "EUR")
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
    i <- 1
    for (c in chr)
    {
        files <- paste0("chr", c, "_", p, "_", pops)
        cat(paste0(c, "..."))
        for (file in files) {
            file_info <- search_for_files(file)
            if (is.null(file_info)) {
                next
            }
            dt_list[[i]] <- fread(cmd = ifelse(file_info$gz,
                paste0("gzcat ", file_info$gene_file),
                paste0("cat ", file_info$gene_file)))
            dt_tmp <- add_N(file_info, dt_list[[i]])
            dt_list[[i]] <- dt_tmp$dt_gene %>% filter(Group != "Cauchy")
            i <- i+1
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
        dt_meta[[test]] <- list()
        Pvalue_col <- ifelse(test == "SKAT-O", "Pvalue", paste0("Pvalue_", test))

        # Weighted Fisher's meta-analysis of p-values
        dt_meta[[test]][["weighted Fisher"]] <- run_weighted_fisher(
            dt %>% group_by(Region, Group, max_MAF),
            "N_eff", Pvalue_col, "Pvalue",
            two_tail = ifelse(test == "Burden", TRUE, FALSE),
            input_beta = ifelse(test == "Burden", "BETA_Burden", NULL)) %>% 
        mutate(Stat = NA, type="weighted Fisher")

        # Stouffer's Z - Make sure P-values match, Stat= weighted_Z_Burden_Stouffer
        dt_meta[[test]][["Stouffer"]] <- run_stouffer(
            dt %>% group_by(Region, Group, max_MAF),
            "N_eff", "Stat", Pvalue_col, "Pvalue",
            two_tail = ifelse(test == "Burden", TRUE, FALSE),
            input_beta = ifelse(test == "Burden", "BETA_Burden", NULL)) %>% 
            mutate(type="Stouffer")
        
        # Cauchy combination meta-analysis, Stat = CCT_Burden
        dt_meta[[test]][["Cauchy"]] <- run_cauchy(
            dt %>% group_by(Region, Group, max_MAF),
            "N_eff", "Stat", Pvalue_col, "Pvalue") %>% 
        select(-number_of_pvals) %>% mutate(type="Cauchy")

        dt_meta[[test]] <- rbindlist(dt_meta[[test]], use.names=TRUE) %>% mutate(class=test)
    }
    dt_meta <- rbindlist(dt_meta)
    p_trim <- gsub("_$", "", str_trim(gsub("[[:space:]_]+", "\\_", p)))
    fwrite(dt_meta, file=paste0("meta_analysis/", p_trim, "_meta_analysis_results_testing.tsv.gz"), sep='\t')
}
