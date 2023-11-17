#!/usr/bin/env Rscript
library(data.table)
library(argparse)
library(dplyr)

CSQ_CODING_HIGH_IMPACT <- c(
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained", 
    "frameshift_variant", "stop_lost")

CSQ_CODING_MEDIUM_IMPACT <- c(
    "start_lost", "initiator_codon_variant", "transcript_amplification", "inframe_insertion",
    "inframe_deletion", "missense_variant", "protein_altering_variant", "splice_region_variant")

CSQ_CODING_LOW_IMPACT <- c(
    "incomplete_terminal_codon_variant", "start_retained_variant", "stop_retained_variant",
    "synonymous_variant", "coding_sequence_variant")

CSQ_NON_CODING <- c(
    "mature_miRNA_variant", "5_prime_UTR_variant", "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant", "non_coding_exon_variant", "intron_variant",
    "NMD_transcript_variant", "non_coding_transcript_variant", "nc_transcript_variant",
    "upstream_gene_variant", "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification",
    "TF_binding_site_variant", "regulatory_region_ablation", "regulatory_region_amplification",
    "feature_elongation", "regulatory_region_variant", "feature_truncation", "intergenic_variant")

# DEV: include the spliceAI annotation information in bins

main <- function(args)
{
    most_deleterious <- function(annotation) {
        case_when(
          "pLoF" %in% annotation ~ "pLoF",
          "damaging_missense_or_protein_altering" %in% annotation ~ "damaging_missense_or_protein_altering",
          "other_missense_or_protein_altering" %in% annotation ~ "other_missense_or_protein_altering",
          "synonymous" %in% annotation ~ "synonymous",
          .default = "non_coding"
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

    # Ensure that all variants above MAX_AF 1% in gnomAD have been removed - do it with and without this.
    AC_path <- args$AC_path
    # python brava_create_annot.py --vep ukb_wes_450k.qced.chr21.vep_processed.txt --spliceai ukb_wes_450k.qced.v6.sites_only.21.all.vcf --out_file test
    vep_spliceAI_path <- args$vep_spliceAI_processed
    
    dt_brava_annot <- fread(vep_spliceAI_path, key=c("SNP_ID"))
    dt_AC <- fread(AC_path)
    dt_AC[, SNP_ID := gsub("chr", "", SNP)]
    setkey(dt_AC, "SNP_ID")
    dt_AC[, AC_A1 := 2*`C(HOM A1)` + `C(HET)`]
    dt_AC[, AC_A2 := 2*`C(HOM A2)` + `C(HET)`]
    dt_AC[, MAC := pmin(AC_A1, AC_A2)]
    dt_AC[, check := `C(HOM A1)` + `C(HOM A2)` + `C(HET)` + `C(MISSING)`]
    dt_AC <- dt_AC[MAC > 0]
    n_samples <- dt_AC$check[1]
    dt_AC[, MAF := MAC/(2*n_samples)]
    dt_AC[, MAF_bin := cut(
        MAF,
        breaks=c(0,(1/n_samples + 1/(2*n_samples)), 0.001, 0.01, 0.05, 0.5, 1),
        labels=c("singletons", "<0.01%", "0.01-0.1%", "1-5%", ">5%", ">50%"),
        include.lowest=TRUE)
    ]
    dt_variant_gene_AC <- merge(dt_AC, dt_brava_annot)
    
    # Transcript level - include variants multiple times where they appear on different transcripts
    dt_summary <- dt_variant_gene_AC %>% group_by(annotation, MAF_bin) %>% summarise(average_AC = sum(MAC))
    dt_summary <- dt_summary %>% mutate(average_count = average_AC/n_samples)
    fwrite(dt_summary, file=paste0(args$out, ".BRaVa_annotations_transcript_summary.tsv.gz"), sep="\t", quote=FALSE)

    dt_variant_AC <- dt_variant_gene_AC %>% group_by(SNP_ID) %>% summarise(annotation = most_deleterious(annotation))
    dt_variant_AC <- data.table(dt_variant_AC, key="SNP_ID")
    dt_variant_AC <- merge(dt_variant_AC, unique(dt_variant_gene_AC %>% select(c("SNP_ID","MAC", "MAF_bin"))))
    
    # Variant level - just include the most deleterious annotation for a variant
    dt_summary <- dt_variant_AC %>% group_by(annotation, MAF_bin) %>% summarise(average_AC = sum(MAC))
    dt_summary <- dt_summary %>% mutate(average_count = average_AC/n_samples)
    fwrite(dt_summary, file=paste0(args$out, ".BRaVa_annotations_variant_summary.tsv.gz"), sep="\t", quote=FALSE)

    # Create summary tables of the result and send to disk.
    # Now do the same thing, but more fine grain, asking what is the distribution of counts of classes of variation output by VEP
    # Use the gnomAD definitions for ordering of deleteriousness to define the class of variant that is annotated.

    names_parse <- grep("&", names(table(dt_variant_gene_AC$CSQ)), value=TRUE)
    mapping <- list()
    names_no_parse <- setdiff(names(table(dt_variant_gene_AC$CSQ)), names_parse)
    mapping[names_no_parse] <- names_no_parse
    for (name in names_parse) {
        CSQs <- strsplit(name, split="&")[[1]]
        mapping[[name]] <- case_when(
            any(CSQs %in% CSQ_CODING_HIGH_IMPACT) ~ ifelse(
                sum(CSQs %in% CSQ_CODING_HIGH_IMPACT) == 1,
                CSQs[CSQs %in% CSQ_CODING_HIGH_IMPACT], 
                paste(CSQs[CSQs %in% CSQ_CODING_HIGH_IMPACT], collapse="&")),
            any(CSQs %in% CSQ_CODING_MEDIUM_IMPACT) ~ ifelse(
                sum(CSQs %in% CSQ_CODING_MEDIUM_IMPACT) == 1,
                CSQs[CSQs %in% CSQ_CODING_MEDIUM_IMPACT],
                paste(CSQs[CSQs %in% CSQ_CODING_MEDIUM_IMPACT], collapse="&")),
            any(CSQs %in% CSQ_CODING_LOW_IMPACT) ~ ifelse(
                sum(CSQs %in% CSQ_CODING_LOW_IMPACT) == 1,
                CSQs[CSQs %in% CSQ_CODING_LOW_IMPACT],
                paste(CSQs[CSQs %in% CSQ_CODING_LOW_IMPACT], collapse="&")),
            any(CSQs %in% CSQ_NON_CODING) ~ ifelse(
                sum(CSQs %in% CSQ_NON_CODING) == 1,
                CSQs[CSQs %in% CSQ_NON_CODING],
                paste(CSQs[CSQs %in% CSQ_NON_CODING], collapse="&")),
            .default = NA
        )
    }

    # Place these within the CSQs column
    dt_variant_gene_AC[, CSQs := unlist(mapping[CSQ])]

    # For the remainders, randomly assign them
    set.seed(42)
    dt_variant_gene_AC[, CSQs := ifelse(grepl("&", CSQs), sample(strsplit(CSQs, split="&")[[1]], 1), CSQs)]

    dt_summary <- dt_variant_gene_AC %>% group_by(CSQs, MAF_bin) %>% summarise(average_AC = sum(MAC))
    dt_summary <- dt_summary %>% mutate(average_count = average_AC/n_samples)
    fwrite(dt_summary, file=paste0(args$out, ".vep_annotations_transcript_summary.tsv.gz"), sep="\t", quote=FALSE)

    dt_variant_AC <- dt_variant_gene_AC %>% group_by(SNP_ID) %>% summarise(CSQ = most_deleterious_vep(CSQs))

    # Take the result and restrict to a single annotation per variant
    dt_variant_AC <- data.table(dt_variant_AC, key="SNP_ID")
    dt_variant_AC <- merge(dt_variant_AC, unique(dt_variant_gene_AC %>% select(c("SNP_ID","MAC", "MAF_bin"))))
    dt_variant_AC[, CSQs := ifelse(grepl("&", CSQ), sample(strsplit(CSQ, split="&")[[1]], 1), CSQ)]

    # Variant level - just include the most deleterious annotation for a variant
    dt_summary <- dt_variant_AC %>% group_by(CSQs, MAF_bin) %>% summarise(average_AC = sum(MAC))
    dt_summary <- dt_summary %>% mutate(average_count = average_AC/n_samples)
    fwrite(dt_summary, file=paste0(args$out, ".vep_annotations_variant_summary.tsv.gz"), sep="\t", quote=FALSE)
}


# Add arguments
parser <- ArgumentParser()
parser$add_argument("--AC_path", default=NULL, required=TRUE, help="Path to allele count information, output from plink after calling --freqx")
parser$add_argument("--vep_spliceAI_processed", default=NULL, required=TRUE, help="Path to the processed VEP file with spliceAI information included. This is the 'long' output of brava_create_annot.py")
parser$add_argument("--out", default=NULL, required=TRUE, help="Output filepath")
args <- parser$parse_args()

main(args)
