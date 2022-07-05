library(data.table)
library(dplyr)

create_pheno_dt <- function(TRANCHE)
{
    PHENOFILE <- paste0("/well/lindgren/UKBIOBANK/dpalmer/ukb_wes_phenotypes/", TRANCHE, '/QC_phenotypes.tsv.gz')
    
    dt <- fread(cmd=paste('zcat', PHENOFILE)) %>% 
        mutate(s=ID) %>% 
        select(
            s, ukbb_centre, Submitted_Gender, Inferred_Gender, genotyping_array, sequencing_batch,
            self_report_ethnicity, genetic_eur_oct2021, genetic_eur_no_fin_oct2021,
            genetically_european, `in_white_British_ancestry_subset`)
    
    # Rename to obtain the non-Europeans
    dt <- dt %>% mutate(
        white_british = ifelse(`in_white_British_ancestry_subset` == 1, "White-British", "Non white-British"),
        genetic_eur_no_fin_oct2021 = ifelse(genetic_eur_no_fin_oct2021, "NFE", "Non-NFE"),
        genetic_eur_oct2021 = ifelse(genetic_eur_oct2021, "European", "Non-European"),
        sequencing_batch = ifelse(sequencing_batch == 1, "Batch 1", "Batch 2"),
        ukbb_centre = factor(ukbb_centre),
        self_report_ethnicity = factor(self_report_ethnicity)
        )

    # Re-level manually.
    dt$self_report_ethnicity <- recode_factor(
        dt$self_report_ethnicity,
        `1` = "White",
        `1001` = "British",
        `2001` = "White and Black Caribbean",
        `3001` = "Indian",
        `4001` = "Caribbean",
        `2` = " Mixed",
        `1002` = "Irish",
        `2002` = "White and Black African",
        `3002` = "Pakistani",
        `4002` = "African",
        `3` = " Asian or Asian British",
        `1003` = "Any other white background",
        `2003` = "White and Asian",
        `3003` = "Bangladeshi",
        `4003` = "Any other Black background",
        `4` = " Black or Black British",
        `2004` = "Any other mixed background",
        `3004` = "Any other Asian background",
        `5` = "Chinese",
        `6` = "Other ethnic group",
        `-1` = "Do not know",
        `-3` = "Prefer not to answer"
        )

    dt$ukbb_centre <- recode_factor(
        dt$ukbb_centre,
        `11012` =  "Barts",
        `11021` =  "Birmingham",
        `11011` =  "Bristol",
        `11008` =  "Bury",
        `11003` =  "Cardiff",
        `11024` =  "Cheadle (revisit)",
        `11020` =  "Croydon",
        `11005` =  "Edinburgh",
        `11004` =  "Glasgow",
        `11018` =  "Hounslow",
        `11010` =  "Leeds",
        `11016` =  "Liverpool",
        `11001` =  "Manchester",
        `11017` =  "Middlesborough",
        `11009` =  "Newcastle",
        `11013` =  "Nottingham",
        `11002` =  "Oxford",
        `11007` =  "Reading",
        `11014` =  "Sheffield",
        `10003` =  "Stockport (pilot)",
        `11006` =  "Stoke",
        `11022` =  "Swansea",
        `11023` =  "Wrexham",
        `11025` =  "Cheadle (imaging)",
        `11026` =  "Reading (imaging)",
        `11027` =  "Newcastle (imaging)",
        `11028` =  "Bristol (imaging)"
        )
    
    dt <- data.table(dt)
    setkey(dt, 's')
    return(dt)
}