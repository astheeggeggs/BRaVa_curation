library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(latex2exp)

source("~/Repositories/BRaVa_curation/QC/utils/pretty_plotting.r")
data_directory <- "~/Repositories/BRaVa_curation/data"
setwd(data_directory)

ribbon_p <- 0.95
files <- dir(
	"meta_analysis", pattern="*meta_analysis_results.tsv.gz", full.names=TRUE
	)

pdf(file="~/Repositories/BRaVa_curation/plots/meta_analysis/meta_analysis_qq.pdf", width=8, height=4)
for (file in files) {
    phe <- gsub(".*/(.*)meta_analysis_results.tsv.gz", "\\1", file)
    phe_plot <- gsub("_", " ", gsub("_$", "", str_trim(gsub("[[:space:]_]+", "\\_", phe))))
    cat(paste0(phe_plot, "...\n"))
    dt_meta <- fread(cmd = paste("gzcat", file)) %>% mutate(Pvalue = -log10(Pvalue))
    dt_meta_to_plot <- data.table(
        dt_meta %>% 
        arrange(class, type, Group, max_MAF, desc(Pvalue)) %>%
        select(-c(class, type, Group, max_MAF)),
        dt_meta %>% 
        group_by(class, type, Group, max_MAF) %>% 
        arrange(desc(Pvalue)) %>% 
        summarize(
            Pvalue_expected = -log10(seq(1, n())/(n() + 1)),
            clower = -log10(qbeta(p = (1 - ribbon_p) / 2, shape2 = n():1, shape1 = 1:n())),
            cupper = -log10(qbeta(p = (1 + ribbon_p) / 2, shape2 = n():1, shape1 = 1:n()))
    	)
    )

    max_MAF_groups <- setdiff(unique(dt_meta_to_plot$max_MAF), "Cauchy")
    groups <- setdiff(unique(dt_meta_to_plot$Group), "Cauchy")

    for (m in max_MAF_groups) {
        for (g in groups) {
            max_MAF_plot <- ifelse(
                grepl("e", as.character(m)),
                gsub("([0-9\\.]+)e(-)*0*([1-9][0-9]*)", "$\\1\\\\times 10^{\\2\\3}$", as.character(m)),
                paste0("$", as.character(m), "$"))
            variant_class_plot <- gsub("_", " ", gsub("[\\|;]", ", ", c))
            p <- create_pretty_qq_plot(
                plot_title=phe_plot,
                plot_subtitle=TeX(paste0(variant_class_plot, "; max MAF = ", max_MAF_plot)),
                cex_labels=2,
                dt_meta_to_plot %>% filter(Group==g, max_MAF==m),
                aes(x=Pvalue_expected, y=Pvalue, color=class),
                save_figure=FALSE,
                x_label=TeX("$-\\log_{10}(P_{expected})$"), 
                y_label=TeX("$-\\log_{10}(P_{observed})$"),
                key_cols=c("class", "Pvalue"),
                aes_ribbon = aes(ymin=clower, ymax=cupper),
                width=170,
                height=120,
                by_chr=FALSE,
                print_p=FALSE
            ) + facet_wrap(~type)
            print(p)
        }
    }

    p <- create_pretty_qq_plot(
        plot_title=phe_plot,
        plot_subtitle="Cauchy",
        cex_labels=2,
        dt_meta_to_plot %>% filter(Group == "Cauchy"),
        aes(x=Pvalue_expected, y=Pvalue, color=class),
        save_figure=FALSE,
        x_label=TeX("$-\\log_{10}(P_{expected})$"), 
        y_label=TeX("$-\\log_{10}(P_{observed})$"),
        key_cols=c("class", "Pvalue"),
        aes_ribbon = aes(ymin=clower, ymax=cupper),
        width=170,
        height=120,
        by_chr=FALSE,
        print_p=FALSE) + facet_wrap(~type)
    print(p)
}
dev.off()
