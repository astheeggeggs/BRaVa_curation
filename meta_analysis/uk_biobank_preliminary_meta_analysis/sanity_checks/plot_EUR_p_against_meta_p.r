library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(latex2exp)

source("~/Repositories/BRaVa_curation/meta_analysis/meta_analysis_utils.r")

data_directory <- "~/Repositories/BRaVa_curation/data"
setwd(data_directory)

chr <- c(seq(1,22), "X")
pops <- "EUR"
max_MAF_groups <- c(1e-4, 1e-3, 1e-2)
groups <- c(
    "damaging_missense_or_protein_altering", 
    "other_missense_or_protein_altering",
    "pLoF",
    "pLoF;damaging_missense_or_protein_altering",
    "pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous",
    "synonymous"
)
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

pdf(file="~/Repositories/BRaVa_curation/plots/meta_analysis/meta_analysis_EUR_comparison.pdf",
    width=6, height=6)
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
            dt_tmp <- fread(cmd = ifelse(file_info$gz,
                paste0("gzcat ", file_info$gene_file),
                paste0("cat ", file_info$gene_file)))
            dt_list[[file]] <- dt_tmp %>% filter(Group != "Cauchy")
            dt_list_cauchy[[file]] <- dt_tmp %>% filter(Group == "Cauchy")
        }
    }

    dt <- rbindlist(dt_list) %>% mutate(max_MAF = as.character(max_MAF))
    dt_cauchy <- rbindlist(dt_list_cauchy) %>% mutate(max_MAF = as.character(max_MAF))
    p_trim <- gsub("_$", "", str_trim(gsub("[[:space:]_]+", "\\_", p)))
    dt_meta <- fread(cmd=paste0("gzcat meta_analysis/", p_trim, "_meta_analysis_results.tsv.gz"))

    dt <- melt(dt, measure.vars=c("Pvalue_Burden", "Pvalue_SKAT", "Pvalue"),
        variable.name = "class", value.name = "Pvalue_EUR")
    dt <- dt %>% mutate(
        class = ifelse(class == "Pvalue", "SKAT-O",
            ifelse(class == "Pvalue_SKAT", "SKAT", "Burden")
            )
        )
    setkeyv(dt, c("Region", "Group", "max_MAF", "class"))
    setkeyv(dt_meta, c("Region", "Group", "max_MAF", "class"))
    dt_merge <- merge(dt, dt_meta)

    for (m in max_MAF_groups) {
        for (g in groups) {
            pl <- ggplot(dt_merge %>% filter(Group == g, max_MAF == m),
                aes(x=-log10(Pvalue_EUR), y=-log10(Pvalue))) + 
            geom_bin_2d(bins=100) + geom_abline(intercept=0,slope=1,col='red',lwd=0.3) +
            facet_grid(rows=vars(class), cols=vars(type)) +
            theme_bw() +
            labs(
                x=TeX("$-\\log_{10}(P_{EUR})$"),
                y=TeX("$-\\log_{10}(P)$"),
                title=gsub("_", " ", p_trim),
                subtitle=g
                )
            print(pl)
        }
    }

    dt_cauchy <- melt(dt_cauchy, measure.vars=c("Pvalue_Burden", "Pvalue_SKAT", "Pvalue"),
        variable.name = "class", value.name = "Pvalue_EUR")
    dt_cauchy <- dt_cauchy %>% mutate(
        class = ifelse(class == "Pvalue", "SKAT-O",
            ifelse(class == "Pvalue_SKAT", "SKAT", "Burden")
            )
        ) %>% filter(Group == "Cauchy") %>% select(-c("max_MAF", "Group"))
    dt_meta <- dt_meta %>% filter(Group == "Cauchy") %>% select(-c("max_MAF", "Group"))
    setkeyv(dt_cauchy, c("Region", "class"))
    setkeyv(dt_meta, c("Region", "class"))
    dt_merge <- merge(dt_cauchy, dt_meta)

    pl <- ggplot(dt_merge,
        aes(x=-log10(Pvalue_EUR), y=-log10(Pvalue))) + 
    geom_bin_2d(bins=100) + geom_abline(intercept=0,slope=1,col='red',lwd=0.3) +
    facet_grid(rows=vars(class), cols=vars(type)) +
    theme_bw() +
    labs(
        x=TeX("$-\\log_{10}(P_{EUR})$"),
        y=TeX("$-\\log_{10}(P)$"),
        title=gsub("_", " ", p_trim),
        subtitle=g
        )
    print(pl)
}
dev.off()
