library(data.table)
library(dplyr)
library(googlesheets4)

# Vocabulary download - ICD10-CM and SNOMED CT from Athena
# 70	CDM 5	ICD10CM	International Classification of Diseases, Tenth Revision, Clinical Modification (NCHS)
# 1	CDM 5	SNOMED	Systematic Nomenclature of Medicine - Clinical Terms (IHTSDO)

dt_concept_1 <- fread("~/Downloads/vocabulary_download_v5_{27f0a1b2-ffd9-42c0-a7d6-6dfae974b856}_1661530459193/CONCEPT.csv")
dt_relationship <- fread("~/Downloads/vocabulary_download_v5_{27f0a1b2-ffd9-42c0-a7d6-6dfae974b856}_1661530459193/CONCEPT_RELATIONSHIP.csv")

# Merge the concept table and relationship tables together

dt_concept_2 <- dt_concept_1

names(dt_concept_1) <- paste0(names(dt_concept_1), "_1")
names(dt_concept_2) <- paste0(names(dt_concept_2), "_2")

setkey(dt_concept_1, "concept_id_1")
setkey(dt_concept_2, "concept_id_2")

setkey(dt_relationship, "concept_id_1")

dt_relationship <- merge(dt_relationship, dt_concept_1)
setkey(dt_relationship, "concept_id_2")

dt_relationship <- merge(dt_relationship, dt_concept_2)

dt_relationship <- dt_relationship %>% select(
	concept_name_1, domain_id_1, vocabulary_id_1, concept_code_1,
	concept_name_2, domain_id_2, vocabulary_id_2, concept_code_2, 
	relationship_id)

# Remove the '.' from the Athena mapping
dt_relationship[, concept_code_1 := gsub("\\.", "", concept_code_1)]
dt_relationship[, concept_code_2 := gsub("\\.", "", concept_code_2)]

dt_relationship_SNOMED_ICD <- dt_relationship %>% filter(
	((vocabulary_id_1 == "ICD10CM") & (vocabulary_id_2 == "SNOMED"))
	)
# # Remove mappings that are clearly wrong
# dt_relationship_SNOMED_ICD %>% filter(grepl("^[A-Z]", concept_code_2))

dt_relationship_SNOMED_SNOMED <- dt_relationship %>% filter(
	((vocabulary_id_1 == "SNOMED") & (vocabulary_id_2 == "SNOMED"))
	)

# Loop over using the grepping
dt <- read_sheet("https://docs.google.com/spreadsheets/d/1YqdSyxf2OyoIYvLnDVj7NmbpebtppsgyJSq18gkVAWI/edit#gid=1716081249", sheet=4, skip=3)
cols <- c(
	"Description",
	"ICD10_control_exclude grep (assuming hierarchy, no .)",
	"ICD10_case_include grep (assuming hierarchy, no .)")
dt <- dt[, cols, with=FALSE]

dt_SNOMED <- dt
dt_SNOMED$SNOMED_case_include <- NA
dt_SNOMED$SNOMED_control_exclude <- NA

for (row in seq(1, nrow(dt_SNOMED))) {
	print(dt_SNOMED$Description[row])
	regexp_case <- paste0(" ", dt_SNOMED$`ICD10_case_include grep (assuming hierarchy, no .)`[row], " ")
	regexp_control <- paste0(" ", dt_SNOMED$`ICD10_control_exclude grep (assuming hierarchy, no .)`[row], " ")
	
	SNOMED_case_include_codes <- paste0(unique(unlist(dt_relationship_SNOMED_ICD %>% filter(grepl(regexp_case, paste0(" ", dt_relationship_SNOMED_ICD$concept_code_1, " "))) %>% select(concept_code_2))), collapse=" | ")
	dt_SNOMED$SNOMED_case_include[row] <- SNOMED_case_include_codes
	print(dt_SNOMED$SNOMED_case_include[row])

	if (!is.na(regexp_control)) {
		SNOMED_control_exclude_codes <- paste0(unique(unlist(dt_relationship_SNOMED_ICD %>% filter(grepl(regexp_control, paste0(" ", dt_relationship_SNOMED_ICD$concept_code_1, " "))) %>% select(concept_code_2))), collapse=" | ")
		dt_SNOMED$SNOMED_control_exclude[row] <- SNOMED_control_exclude_codes
		print(dt_SNOMED$SNOMED_control_exclude[row])
	}
}

# Now loop through and include all SNOMED codes subsumed by the SNOMED codes available for each phenotype
dt_relationship_SNOMED_SNOMED <- dt_relationship_SNOMED_SNOMED %>% filter(relationship_id %in% c("Maps to", "Subsumes"))
dt_SNOMED$SNOMED_case_include_broad <- NA
dt_SNOMED$SNOMED_control_exclude_broad <- NA

for (row in seq(1, nrow(dt_SNOMED))) {
	print(dt_SNOMED$Description[row])
	regexp_case <- paste0(" ", dt_SNOMED$SNOMED_case_include[row], " ")
	regexp_control <- paste0(" ", dt_SNOMED$SNOMED_control_exclude[row], " ")
	
	SNOMED_case_include_codes <- unique(unlist(dt_relationship_SNOMED_SNOMED %>% filter(grepl(regexp_case, paste0(" ", dt_relationship_SNOMED_SNOMED$concept_code_1, " "))) %>% select(concept_code_2)))
	SNOMED_case_include_codes <- paste0(unique(SNOMED_case_include_codes, strsplit(gsub(" ", "", regexp_case), split="\\|")[[1]]), collapse=" | ")
	dt_SNOMED$SNOMED_case_include_broad[row] <- SNOMED_case_include_codes
	print(dt_SNOMED$SNOMED_case_include_broad[row])

	if (!is.na(regexp_control)) {
		SNOMED_control_exclude_codes <- unique(unlist(dt_relationship_SNOMED_SNOMED %>% filter(grepl(regexp_control, paste0(" ", dt_relationship_SNOMED_SNOMED$concept_code_1, " "))) %>% select(concept_code_2)))
		SNOMED_control_exclude_codes <- paste0(unique(SNOMED_control_exclude_codes, strsplit(gsub(" ", "", regexp_control), split="\\|")[[1]]), collapse=" | ")
		dt_SNOMED$SNOMED_control_exclude_broad[row] <- SNOMED_control_exclude_codes
		print(dt_SNOMED$SNOMED_case_include_broad[row])
	}
}

fwrite(dt_SNOMED, sep="\t", file="data/SNOMED_brava_ICD.tsv")
