#!/usr/bin/env Rscript
library(data.table)
library(argparse)
library(dplyr)

CSQ_CODING_HIGH_IMPACT <- c(
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
    "stop_gained", "frameshift_variant", "stop_lost")

CSQ_CODING_MEDIUM_IMPACT <- c(
    "start_lost", "initiator_codon_variant", "transcript_amplification",
    "inframe_insertion", "inframe_deletion", "missense_variant",
    "protein_altering_variant", "splice_region_variant")

CSQ_CODING_LOW_IMPACT <- c(
    "incomplete_terminal_codon_variant", "start_retained_variant",
    "stop_retained_variant", "synonymous_variant", "coding_sequence_variant")

CSQ_NON_CODING <- c(
    "mature_miRNA_variant", "5_prime_UTR_variant", "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant", "non_coding_exon_variant",
    "intron_variant", "NMD_transcript_variant", "non_coding_transcript_variant",
    "nc_transcript_variant", "upstream_gene_variant", "downstream_gene_variant",
    "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant",
    "regulatory_region_ablation", "regulatory_region_amplification",
    "feature_elongation", "regulatory_region_variant", "feature_truncation",
    "intergenic_variant")

main <- function(args)
{
    most_deleterious <- function(annotation) {
        case_when(
          "pLoF" %in% annotation ~ "pLoF",
          "damaging_missense_or_protein_altering" %in% 
          annotation ~ "damaging_missense_or_protein_altering",
          "other_missense_or_protein_altering" %in% 
          annotation ~ "other_missense_or_protein_altering",
          "synonymous" %in% annotation ~ "synonymous",
          "non_coding" %in% annotation ~ "non_coding",
          .default = NA
        )
    }

    most_deleterious_vep <- function(CSQ) {
        CSQ <- unique(CSQ)
        case_when(
            any(CSQ %in% CSQ_CODING_HIGH_IMPACT) ~ ifelse(
                sum(CSQ %in% CSQ_CODING_HIGH_IMPACT) == 1,
                CSQ[CSQ %in% CSQ_CODING_HIGH_IMPACT], 
                paste(CSQ[CSQ %in% CSQ_CODING_HIGH_IMPACT], collapse="&")),
            any(CSQ %in% CSQ_CODING_MEDIUM_IMPACT) ~ ifelse(
                sum(CSQ %in% CSQ_CODING_MEDIUM_IMPACT) == 1,
                CSQ[CSQ %in% CSQ_CODING_MEDIUM_IMPACT],
                paste(CSQ[CSQ %in% CSQ_CODING_MEDIUM_IMPACT], collapse="&")),
            any(CSQ %in% CSQ_CODING_LOW_IMPACT) ~ ifelse(
                sum(CSQ %in% CSQ_CODING_LOW_IMPACT) == 1,
                CSQ[CSQ %in% CSQ_CODING_LOW_IMPACT],
                paste(CSQ[CSQ %in% CSQ_CODING_LOW_IMPACT], collapse="&")),
            any(CSQ %in% CSQ_NON_CODING) ~ ifelse(
                sum(CSQ %in% CSQ_NON_CODING) == 1,
                CSQ[CSQ %in% CSQ_NON_CODING],
                paste(CSQ[CSQ %in% CSQ_NON_CODING], collapse="&")),
            .default = NA
        )
    }

    create_MAC_and_MAF_summary <- function(
        dt, annotation_group, fout, spliceAI_bins=TRUE)
    {
        annotation_group <- enquo(annotation_group)

        dt_summary_MAF <- dt %>% group_by(!!annotation_group, MAF_bin) %>% 
            summarise(average_AC = sum(MAC), variant_count = n())
        dt_summary_MAC <- dt %>% group_by(!!annotation_group, MAC_bin) %>% 
            summarise(average_AC = sum(MAC), variant_count = n())
        dt_summary_MAF <- dt_summary_MAF %>% rename(bin=MAF_bin) %>% 
            mutate(bin_type="MAF")
        dt_summary_MAC <- dt_summary_MAC %>% rename(bin=MAC_bin) %>% 
            mutate(bin_type="MAC")

        dt_summary_MAF <- data.table(dt_summary_MAF)
        dt_summary_MAC <- data.table(dt_summary_MAC)
        dt_summary <- merge(dt_summary_MAF, dt_summary_MAC, all=TRUE)

        if (spliceAI_bins) {
            dt_summary_spliceAI <- dt %>%
                group_by(!!annotation_group, spliceAI_bin) %>% 
                summarise(average_AC = sum(MAC), variant_count = n())
            dt_summary_spliceAI <- dt_summary_spliceAI %>%
            rename(bin=spliceAI_bin) %>% mutate(bin_type="spliceAI")
            dt_summary <- merge(dt_summary, dt_summary_spliceAI, all=TRUE)
        }

        dt_summary <- dt_summary %>% 
            mutate(average_count = average_AC/n_samples)

        fwrite(dt_summary, file=fout, sep="\t", quote=FALSE)
        cat("created summary split by annotation\n")
        return(dt_summary)
    }

    # Ensure that all variants above MAX_AF 1% in gnomAD have been removed 
    # do it with and without this.
    AC_path <- args$AC_path
    # python brava_create_annot.py 
    #       --vep ukb_wes_450k.qced.chr21.vep_processed.txt 
    #       --spliceai ukb_wes_450k.qced.v6.sites_only.21.all.vcf 
    #       --out_file test
    vep_spliceAI_path <- args$vep_spliceAI_processed
    
    dt_brava_annot <- fread(vep_spliceAI_path, key=c("SNP_ID"))
    if (!all(dt_brava_annot$annotation %in% c(
        "pLoF",
        "damaging_missense_or_protein_altering",
        "other_missense_or_protein_altering",
        "synonymous",
        "non_coding"))) {
        stop("An annotation is present in the file which is not in the following:
            - pLoF
            - damaging_missense_or_protein_altering
            - other_missense_or_protein_altering
            - synonymous
            - non_coding")
    }
    if (all(is.na(dt_brava_annot$max_DS))) {
        stop("splice AI information has not been incorporated.
            All entries of the max_DS column are empty. 
            Check variant ID merging between spliceAI and VEP information. 
            Perhaps chrCHR:POS:REF:ALT vs CHR:POS:REF:ALT.")
    }
    dt_brava_annot[, SNP_ID := gsub("chr", "", SNP_ID)]
    dt_AC <- fread(AC_path)

    dt_AC[, SNP_ID := gsub("chr", "", SNP)]
    setkey(dt_AC, "SNP_ID")

    # Check to ensure that all SNP_IDs are off the form CHR:POS:REF:ALT
    AC_correct_format <- all(
        grepl("^[1-9,X]{1}[0-9]*:[0-9]+:[A,C,G,T]+:[A,C,G,T]+", dt_AC$SNP_ID))
    brava_annot_correct_format <- all(
        grepl("^[1-9,X]{1}[0-9]*:[0-9]+:[A,C,G,T]+:[A,C,G,T]+", dt_brava_annot$SNP_ID)
        )
    if (!AC_correct_format | !brava_annot_correct_format) {
        stop(paste("One or both SNP/SNP_ID columns in", AC_path,
            "and", vep_spliceAI_path, "is not formatted to CHR:POS:REF:ALT.",
            " Note that chrs 1-9 cannot be labelled 01, 02, etc"))
    }

    dt_AC[, AC_A1 := 2*`C(HOM A1)` + `C(HET)`]
    dt_AC[, AC_A2 := 2*`C(HOM A2)` + `C(HET)`]
    dt_AC[, MAC := pmin(AC_A1, AC_A2)]
    dt_AC[, check := `C(HOM A1)` + `C(HOM A2)` + `C(HET)` + `C(MISSING)`]
    dt_AC <- dt_AC[MAC > 0]
    n_samples <- dt_AC$check[1]
    # dt_AC[, MAF := MAC/(2*n_samples)]
    dt_AC[, MAF := MAC/(AC_A1 + AC_A2)]
    dt_AC[, MAF_bin := cut(
        MAF,
        breaks=c(0, 0.001, 0.01, 0.05, 0.5, 1),
        labels=c("<0.1%", "0.1-1%", "1-5%", ">5%", ">50%"),
        include.lowest=TRUE)
    ]
    dt_AC[, MAC_bin := cut(
        MAC,
        breaks=c(0, 1, 5, 10, 100, 1000, 10000, Inf),
        labels=c("singletons", "(1,5]", "(5,10]", "(10,100]", "(100,1,000]",
            "(1,000,10,000]", ">10,000"))]
    cols <- c("SNP_ID","MAC", "MAF_bin", "MAC_bin")
    dt_variant_gene_AC <- merge(dt_AC, dt_brava_annot)

    if (args$spliceAI_bins) {
        dt_variant_gene_AC[, spliceAI_bin := cut(
            max_DS,
            breaks=c(0, 0.2, 0.5, 0.8, Inf),
            labels=c("<0.2", "[0.2,0.5)", "[0.5,0.8)", ">0.8"))]
        cols <- c("SNP_ID","MAC", "MAF_bin", "MAC_bin", "spliceAI_bin")
    }
    
    # Transcript level
    # Include variants multiple times where they appear on different transcripts
    dt_summary <- create_MAC_and_MAF_summary(
        dt_variant_gene_AC, annotation,
        paste0(args$out, ".BRaVa_annotations_transcript_summary.tsv.gz"),
        spliceAI_bins=args$spliceAI_bins)

    dt_variant_AC <- dt_variant_gene_AC %>% group_by(SNP_ID) %>% 
        summarise(annotation = most_deleterious(annotation))
    dt_variant_AC <- data.table(dt_variant_AC, key="SNP_ID")
    dt_variant_AC <- merge(dt_variant_AC,
        unique(dt_variant_gene_AC %>% select(all_of(cols))))
    
    # Variant level
    # Just include the most deleterious annotation for a variant
    create_MAC_and_MAF_summary(
        dt_variant_AC, annotation,
        paste0(args$out, ".BRaVa_annotations_variant_summary.tsv.gz"),
        spliceAI_bins=args$spliceAI_bins)

    names_parse <- grep("&", names(table(dt_variant_gene_AC$CSQ)), value=TRUE)
    mapping <- list()
    names_no_parse <- setdiff(names(table(dt_variant_gene_AC$CSQ)), names_parse)
    mapping[names_no_parse] <- names_no_parse
    for (name in names_parse) {
        CSQs <- strsplit(name, split="&")[[1]]
        mapping[[name]] <- most_deleterious_vep(CSQs)
    }

    # Place these within the CSQs column
    dt_variant_gene_AC[, CSQs := unlist(mapping[CSQ])]

    # For the remainder, randomly assign them
    set.seed(42)
    dt_variant_gene_AC[, CSQs := ifelse(
        grepl("&", CSQs), sample(strsplit(CSQs, split="&")[[1]], 1), CSQs)]
    
    # Transcript level
    # include variants multiple times where they appear on different transcripts
    create_MAC_and_MAF_summary(
        dt_variant_gene_AC, CSQs,
        paste0(args$out, ".vep_annotations_transcript_summary.tsv.gz"),
        spliceAI_bins=args$spliceAI_bins)

    # Take the result and restrict to a single annotation per variant
    dt_variant_AC <- dt_variant_gene_AC %>% group_by(SNP_ID) %>% 
        summarise(CSQ = most_deleterious_vep(CSQs))
    dt_variant_AC <- data.table(dt_variant_AC, key="SNP_ID")
    dt_variant_AC <- merge(dt_variant_AC,
        unique(dt_variant_gene_AC %>% select(all_of(cols))))
    dt_variant_AC[, CSQs := ifelse(
        grepl("&", CSQ), sample(strsplit(CSQ, split="&")[[1]], 1), CSQ)]

    # Variant level
    # just include the most deleterious annotation for a variant
    dt_summary <- create_MAC_and_MAF_summary(
        dt_variant_AC, CSQs,
        paste0(args$out, ".vep_annotations_variant_summary.tsv.gz"),
        spliceAI_bins=args$spliceAI_bins)
}

# Add arguments
parser <- ArgumentParser()
parser$add_argument("--AC_path", default=NULL, required=TRUE,
    help="Path to allele count information, output from plink --freqx")
parser$add_argument("--vep_spliceAI_processed", default=NULL, required=TRUE,
    help=paste0("Path to the processed VEP file with spliceAI information ",
        "included. This is the 'long' output of brava_create_annot.py"))
parser$add_argument("--out", default=NULL, required=TRUE,
    help="Output filepath")
parser$add_argument("--spliceAI_bins", default=FALSE, action="store_true")
args <- parser$parse_args()

main(args)
