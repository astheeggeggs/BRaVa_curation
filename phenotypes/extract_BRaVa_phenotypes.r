# install.packages("googlesheets4")
library(data.table)
library(dplyr)
library(ggplot2)
library(googlesheets4)

munge_BRaVa_ICD_proposals <- function() {
	# Download and extract the case and control ICD9/10 codes from the BRaVa nominate phenotypes file
	dt <- read_sheet("https://docs.google.com/spreadsheets/d/1YqdSyxf2OyoIYvLnDVj7NmbpebtppsgyJSq18gkVAWI/edit#gid=1716081249", sheet=4, skip=3)
	cols <- c(
		"Description",
		"ICD10_control_exclude", "ICD10_case_include",
		"ICD9_control_exclude", "ICD9_case_include")
	dt <- dt[, cols, with=FALSE]
	names(dt) <- c("phenotype", "ICD10_control_exclude", "ICD10_case_include", "ICD9_control_exclude", "ICD9_case_include")
	dt <- data.table(dt)

	# Remove any extra spaces
	dt[, ICD10_control_exclude := gsub("[[:space:]]", "", ICD10_control_exclude)]
	dt[, ICD10_case_include := gsub("[[:space:]]", "", ICD10_case_include)]
	dt[, ICD9_control_exclude := gsub("[[:space:]]", "", ICD9_control_exclude)]
	dt[, ICD9_case_include := gsub("[[:space:]]", "", ICD9_case_include)]
	
	# Replace '.'s with '' (it's an equivalent encoding for ICD codes).
	dt[, ICD10_control_exclude := gsub("\\.", "", ICD10_control_exclude)]
	dt[, ICD10_case_include := gsub("\\.", "", ICD10_case_include)]
	dt[, ICD9_control_exclude := gsub("\\.", "", ICD9_control_exclude)]
	dt[, ICD9_case_include := gsub("\\.", "", ICD9_case_include)]

	# Replace commas for pipes to get the regular expression ready.
	dt[, ICD10_control_exclude := gsub(",|;", " | ", ICD10_control_exclude)]
	dt[, ICD10_case_include := gsub(",|;", " | ", ICD10_case_include)]
	dt[, ICD9_control_exclude := gsub(",|;", " | ", ICD9_control_exclude)]
	dt[, ICD9_case_include := gsub(",|;", " | ", ICD9_case_include)]

	# Ensure that the first and final character of the regular expression is a space.
	dt[, ICD10_control_exclude := gsub("^(.*)$", " \\1 ", ICD10_control_exclude)]
	dt[, ICD10_case_include := gsub("^(.*)$", " \\1 ", ICD10_case_include)]
	dt[, ICD9_control_exclude := gsub("^(.*)$", " \\1 ", ICD9_control_exclude)]
	dt[, ICD9_case_include := gsub("^(.*)$", " \\1 ", ICD9_case_include)]

	# Remove asterisks
	dt[, ICD10_control_exclude := gsub("\\*", "", ICD10_control_exclude)]
	dt[, ICD10_case_include := gsub("\\*", "", ICD10_case_include)]
	dt[, ICD9_control_exclude := gsub("\\*", "", ICD9_control_exclude)]
	dt[, ICD9_case_include := gsub("\\*", "", ICD9_case_include)]

	dt$ICD10_control_exclude[which(is.na(dt$ICD10_control_exclude))] <- ""
	dt$ICD10_case_include[which(is.na(dt$ICD10_case_include))] <- ""
	dt$ICD9_control_exclude[which(is.na(dt$ICD9_control_exclude))] <- ""
	dt$ICD9_case_include[which(is.na(dt$ICD9_case_include))] <- ""

	return(dt)
}

extract_ID_and_ICD_UKB <- function(
	phenotype_file = "/well/lindgren-ukbb/projects/ukbb-11867/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844_ukb50009_updateddiagnoses_14012022.csv"
)
{
	get_cols <- function(codes, dt, na.filter=FALSE)
	{
		cols <- c()
		for (code in codes) { cols <- c(cols, grep(paste0("^", code, "\\-"), names(dt), value=TRUE)) }
		return(cols)
	}

	dt <- fread(phenotype_file, na.strings=NULL, nrow=1)

	# Extract the relevant columns for each of the encodings
	ICD10s <- c("41202", "41204", "40006", "40001", "40002")
	ICD9s <- c("41203", "41205", "40013")

	cols <-  get_cols(c(ICD9s, ICD10s), dt)
	select_cols <- rep("character", (length(cols) + 1))
	names(select_cols) <- c("eid", cols)

	# Read in the entire file ensuring these columns are encoded as characters to avoid NA weirdness.
	dt <- fread(phenotype_file, na.strings=NULL, select=select_cols)
	
	# Extract the ICD10 columns
	ICD10_cols <-  get_cols(ICD10s, dt)

	# Combine vector into a single string
	dt[, ICD10_string := do.call(paste,.SD), .SDcols=ICD10_cols]
	# Remove NAs
	dt[, ICD10_string := gsub("NA", "", ICD10_string)]
	# Pad start and end with space
	dt[, ICD10_string := gsub("^(.*)$", " \\1 ", ICD10_string)]
	# Remove trailing whitespace
	dt[, ICD10_string := gsub("( )+", " ", ICD10_string)]
	# Remove the '.'s from the file (there aren't any, but this will be useful for collaborators).
	dt[, ICD10_string := gsub("\\.", "", ICD10_string)]

	# Extract the ICD9 columns
	ICD9_cols <-  get_cols(ICD9s, dt)

	# Combine vector into a single string
	dt[, ICD9_string := do.call(paste,.SD), .SDcols=ICD9_cols]
	# Remove NAs
	dt[, ICD9_string := gsub("NA", "", ICD9_string)]
	# Pad start and end with space
	dt[, ICD9_string := gsub("^(.*)$", " \\1 ", ICD9_string)]
	# Remove trailing whitespace
	dt[, ICD9_string := gsub("( )+", " ", ICD9_string)]

	# Remove the '.'s from the file (there aren't any, but this will be useful for collaborators).
	dt[, ICD9_string := gsub("\\.", "", ICD9_string)]

	if (length(which(is.na(dt$ICD10_string))) > 0) {
		dt$ICD10_string[which(is.na(dt$ICD10_string))] <- ""
	}

	if (length(which(is.na(dt$ICD9_string))) > 0) {
		dt$ICD9_string[which(is.na(dt$ICD9_string))] <- ""
	}

	cols <- c("eid", "ICD9_string", "ICD10_string")
	dt <- dt[, cols, with=FALSE]
	return(dt)
}

extract_case_status <- function(dt_data, dt_query, assume_tree=TRUE) {
	
	cols <- dt_query$phenotype
	dt_data[, (cols) := 0]

	simplify_query <- function(query) {
		# Simplify to speed things up
		query_tmp <- gsub("^ (.*) ", "\\1", query)
		query_tmp <- strsplit(query_tmp, split=" \\| ")[[1]]
		query_tmp <- query_tmp[order(nchar(query_tmp))]
		query <- c()
		while (length(query_tmp) > 0) {
			query <- c(query, query_tmp[1])
			match <- grep(query[length(query)], query_tmp)
			query_tmp <- query_tmp[-match]

		}
		query <- paste0(query, collapse=" | ")
		return(query)
	}

	for (i in 1:length(cols)) {
		# Add a new column to dt_data
		phenotype <- cols[i]
		cat(paste0(phenotype, "...\n"))

		if (dt_query$ICD10_control_exclude[i] != "")
		{
			query <- dt_query$ICD10_control_exclude[i]
			if (assume_tree) {
				query <- simplify_query(query)
				query <- gsub(" \\| ", ".* \\| ", query)
				if (query != "") {
					query <- gsub("$", ".*", query)
				}
				print(query)
			}
			ICD10_where <- grep(query, dt_data$ICD10_string)
			cat(paste0("ICD10 control exclusions: ", length(ICD10_where), "...\n"))
			dt_data[[phenotype]][ICD10_where] <- NA
			dt_query$ICD10_control_exclude[i] <- query
		}

		if (dt_query$ICD9_control_exclude[i] != "")
		{
			query <- dt_query$ICD9_control_exclude[i]
			if (assume_tree) {
				query <- simplify_query(query)
				query <- gsub(" \\| ", ".* \\| ", query)
				if (query != "") {
					query <- gsub("$", ".*", query)
				}
				print(query)
			}
			ICD9_where <- grep(query, dt_data$ICD9_string)
			cat(paste0("ICD9 control exclusions: ", length(ICD9_where), "...\n"))
			dt_data[[phenotype]][ICD9_where] <- NA
			dt_query$ICD9_control_exclude[i] <- query
		}

		cat(paste0("total control exclusions count prior to adding cases: ", sum(is.na(dt_data[[phenotype]])), "\n"))

		if (dt_query$ICD10_case_include[i] != "")
		{
			query <- dt_query$ICD10_case_include[i]
			if (assume_tree) {
				query <- simplify_query(query)
				query <- gsub(" \\| ", ".* \\| ", query)
				if (query != "") {
					query <- gsub("$", ".*", query)
				}
				print(query)
			}
			ICD10_where <- grep(query, dt_data$ICD10_string)
			cat(paste0("ICD10 cases: ", length(ICD10_where), "...\n"))
			dt_data[[phenotype]][ICD10_where] <- 1
			dt_query$ICD10_case_include[i] <- query
		}

		if (dt_query$ICD9_case_include[i] != "")
		{
			query <- dt_query$ICD9_case_include[i]
			if (assume_tree) {
				query <- simplify_query(query)
				query <- gsub(" \\| ", ".* \\| ", query)
				if (query != "") {
					query <- gsub("$", ".*", query)
				}
				print(query)
			}
			ICD9_where <- grep(query, dt_data$ICD9_string)
			cat(paste0("ICD9 cases: ", length(ICD9_where), "...\n"))
			dt_data[[phenotype]][ICD9_where] <- 1
			dt_query$ICD9_case_include[i] <- query
		}

		cat(paste0("total case count: ", sum(dt_data[[phenotype]], na.rm=TRUE), "\n"))
		cat(paste0("final case count: ", sum(dt_data[[phenotype]], na.rm=TRUE), "\n"))
		cat(paste0("final control exclusions count after adding cases: ", sum(is.na(dt_data[[phenotype]])), "\n"))
		cat(paste0("final control count: ", sum(dt_data[[phenotype]] == 0, na.rm=TRUE), "\n"))
	}
	return(list(dt_data=dt_data, dt_query=dt_query))
}

extract_continuous_trait_counts <- function(
	phenotype_file = "/well/lindgren-ukbb/projects/ukbb-11867/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844_ukb50009_updateddiagnoses_14012022.csv",
	biomarker_file = "/well/lindgren/UKBIOBANK/DATA/Biomarker_data/ukb27722.csv",
	manual_curation_mapping = "data/BRaVa_continuous_trait_manual_mappings.tsv",
	superpopulation_labels = "/well/lindgren/UKBIOBANK/dpalmer/superpopulation_assignments/superpopulation_labels.tsv")
{
	get_cols <- function(codes, dt, na.filter=FALSE)
	{
		cols <- c()
		for (code in codes) { cols <- c(cols, grep(paste0("^", code, "\\-"), names(dt), value=TRUE)) }
		return(cols)
	}

	dt_pheno <- fread(phenotype_file, na.strings=NULL, nrow=1)

	# First the columns in the phenotype file
	dt_manual <- fread(manual_curation_mapping) %>% mutate(UKB_code = as.character(UKB_code))
	cols <- (dt_manual %>% filter(UKB_code != ""))$UKB_code
	cols <-  get_cols(cols, dt_pheno)
	pheno_cols <- grep("*\\-0\\.0", cols, value=TRUE)
	select_cols <- c("character", rep("numeric", length(pheno_cols)))
	names(select_cols) <- c("eid", pheno_cols)

	# Read in the entire file ensuring these columns are encoded as characters to avoid NA weirdness.
	dt_pheno <- fread(phenotype_file, select=select_cols)
	setkey(dt_pheno, "eid")

	# Next, the columns in the biomarker file
	dt_biomarker <- fread(biomarker_file, na.strings=NULL, nrow=1)

	cols <- (dt_manual %>% filter(UKB_code != ""))$UKB_code
	cols <-  get_cols(cols, dt_biomarker)
	biomarker_cols <- grep("*\\-0\\.0", cols, value=TRUE)
	select_cols <- c("character", rep("numeric", length(biomarker_cols)))
	names(select_cols) <- c("eid", biomarker_cols)

	# Read in the entire file ensuring these columns are encoded as characters to avoid NA weirdness.
	dt_biomarker <- fread(biomarker_file, select=select_cols)

	# Get population specific counts

	# Merge with 1000G labels - this file is created using 05_estimate_superpopulation.r in the QC folder.
	dt_classify <- fread("/well/lindgren/UKBIOBANK/dpalmer/superpopulation_assignments/superpopulation_labels.tsv")
	dt_classify[, eid:=sample.ID]
	dt_classify[, sample.ID:=NULL]
	setkey(dt_classify, "eid")

	dt_cts_classified <- merge(merge(dt_pheno, dt_classify), dt_biomarker)

	# Finally, take the sum of the BMI, waist, and hip circumference to determine the non-NA sample size for WHR adjusted for BMI.
	dt_cts_classified <- dt_cts_classified %>% mutate(WHR_adjusted_BMI = `48-0.0` + `49-0.0` + `21001-0.0`)

	sum_not_is.na <- function(col) { sum(!is.na(col)) }

	# Split by 1000G label and count the number of non-NA entries
	dt_counts <- dt_cts_classified %>% group_by(classification_strict) %>% summarise(across(c(pheno_cols, biomarker_cols, "WHR_adjusted_BMI"), sum_not_is.na))
	dt_counts_t <- data.table::transpose(dt_counts, keep.names="phenotype", make.names="classification_strict") %>% mutate(phenotype = gsub("-.*", "", phenotype)) %>% rename(UKB_code = phenotype)
	dt_counts_t <- merge(dt_manual, dt_counts_t, all.y=TRUE)

	fwrite(dt_counts_t, file="data/output/UKBB_cts_non_missing_counts.tsv", sep="\t")

	dt_counts <- dt_cts_classified %>% summarise(across(c(pheno_cols, biomarker_cols, "WHR_adjusted_BMI"), sum_not_is.na))
	dt_counts_t <- data.table::transpose(dt_counts, keep.names="phenotype") %>% mutate(phenotype = gsub("-.*", "", phenotype)) %>% rename(UKB_code = phenotype)
	dt_counts_t <- merge(dt_manual, dt_counts_t, all.y=TRUE)

	fwrite(dt_counts_t, file="data/output/UKBB_cts_non_missing_counts_total.tsv", sep="\t")

}

munge_BRaVa_OPCS4_proposals <- function() {
	# Download and extract the case and control ICD9/10 codes from the BRaVa nominate phenotypes file
	dt <- read_sheet("https://docs.google.com/spreadsheets/d/1YqdSyxf2OyoIYvLnDVj7NmbpebtppsgyJSq18gkVAWI/edit#gid=1716081249", sheet=8)
	cols <- c("Procedures", "OPCS codes")

	dt <- dt[, cols, with=FALSE]
	names(dt) <- c("phenotype", "OPCS4_code")
	dt <- data.table(dt)

	# Remove any extra spaces
	dt[, OPCS4_code := gsub("[[:space:]]", "", OPCS4_code)]
	# Replace '.'s with '' (it's an equivalent encoding for OPCS4 codes).
	dt[, OPCS4_code := gsub("\\.", "", OPCS4_code)]
	# Replace commas and semi-colons for pipes to get the regular expression ready.
	dt[, OPCS4_code := gsub(",|;", " | ", OPCS4_code)]
	# Ensure that the first and final character of the regular expression is a space.
	dt[, OPCS4_code := gsub("^(.*)$", " \\1 ", OPCS4_code)]
	# Remove asterisks
	dt[, OPCS4_code := gsub("\\*", "", OPCS4_code)]
	dt$OPCS4_code[which(is.na(dt$IOPCS4_code))] <- ""

	return(dt)
}

extract_ID_and_OPCS4_UKB <- function(
	phenotype_file = "/well/lindgren-ukbb/projects/ukbb-11867/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844_ukb50009_updateddiagnoses_14012022.csv",
	superpopulation_labels = "/well/lindgren/UKBIOBANK/dpalmer/superpopulation_assignments/superpopulation_labels.tsv")
{

	get_cols <- function(codes, dt, na.filter=FALSE)
	{
		cols <- c()
		for (code in codes) { cols <- c(cols, grep(paste0("^", code, "\\-"), names(dt), value=TRUE)) }
		return(cols)
	}

	dt <- fread(phenotype_file, na.strings=NULL, nrow=1)

	# Extract the relevant columns for each of the encodings
	OPCS4s <- c("41200", "41210", "41272")
	cols <-  get_cols(OPCS4s, dt)
	select_cols <- rep("character", (length(cols) + 1))
	names(select_cols) <- c("eid", cols)

	# Read in the entire file ensuring these columns are encoded as characters to avoid NA weirdness.
	dt <- fread(phenotype_file, na.strings=NULL, select=select_cols)

	# Combine vector into a single string
	dt[, OPCS4_string := do.call(paste,.SD), .SDcols=cols]
	# Remove NAs
	dt[, OPCS4_string := gsub("NA", "", OPCS4_string)]
	# Pad start and end with space
	dt[, OPCS4_string := gsub("^(.*)$", " \\1 ", OPCS4_string)]
	# Remove trailing whitespace
	dt[, OPCS4_string := gsub("( )+", " ", OPCS4_string)]
	# Remove the '.'s from the file (there aren't any, but this will be useful for collaborators).
	dt[, OPCS4_string := gsub("\\.", "", OPCS4_string)]

	if (length(which(is.na(dt$OPCS4_string))) > 0) {
		dt$OPCS4_string[which(is.na(dt$OPCS4_string))] <- ""
	}

	cols <- c("eid", "OPCS4_string")
	dt <- dt[, cols, with=FALSE]
	return(dt)
}

extract_procedure_case_status <- function(dt_data, dt_query)
{	
	cols <- dt_query$phenotype
	dt_data[, (cols) := 0]
	
	for (i in 1:length(cols)) {
		# Add a new column to dt_data
		phenotype <- cols[i]
		cat(paste0(phenotype, "...\n"))

		if (dt_query$OPCS4_code[i] != "")
		{
			query <- dt_query$OPCS4_code[i]
			OPCS4_where <- grep(query, dt_data$OPCS4_string)
			cat(paste0("OPCS4 code cases: ", length(OPCS4_where), "...\n"))
			dt_data[[phenotype]][OPCS4_where] <- 1
		}
		cat(paste0("total case count: ", sum(dt_data[[phenotype]], na.rm=TRUE), "\n"))
	}
	return(dt_data)
}

# Extract continous traits
extract_continuous_trait_counts()

dt_query_ICD <- munge_BRaVa_ICD_proposals()
fwrite(dt_query_ICD, file='data/BRaVa_ICD_proposals.tsv', sep='\t', quote=TRUE)
dt_query_ICD <- fread('data/BRaVa_ICD_proposals.tsv')
dt_data <- extract_ID_and_ICD_UKB()
dt_binary_list <- extract_case_status(dt_data, dt_query_ICD)
dt_binary_ICD <- dt_binary_list$dt_data
dt_query_ICD <- dt_binary_list$dt_query
fwrite(dt_query_ICD, file='data/BRaVa_ICD_proposals_grep.tsv', sep='\t', quote=TRUE)
setkey(dt_binary_ICD, "eid")

# Merge with 1000G labels - this file is created using 05_estimate_superpopulation.r in the QC folder.
dt_classify <- fread("/well/lindgren/UKBIOBANK/dpalmer/superpopulation_assignments/superpopulation_labels.tsv")
dt_classify[, eid:=sample.ID]
dt_classify[, sample.ID:=NULL]
setkey(dt_classify, "eid")

dt_binary_classified <- merge(dt_binary_ICD, dt_classify)

dt_query_OPCS4 <- munge_BRaVa_OPCS4_proposals()
fwrite(dt_query_OPCS4, file='data/BRaVa_OPCS4_proposals.tsv', sep='\t', quote=TRUE)
dt_query_OPCS4 <- fread('data/BRaVa_OPCS4_proposals.tsv')
dt_data <- extract_ID_and_OPCS4_UKB()

dt_binary_procedures <- extract_procedure_case_status(dt_data, dt_query_OPCS4)
ICD_and_procedures <- intersect(names(dt_binary_procedures)[-1], names(dt_binary_ICD)[-1])
if (length(ICD_and_procedures) >= 1) {
	where <- which(names(dt_binary_procedures) %in% ICD_and_procedures)
	names(dt_binary_procedures)[where] <- paste0(names(dt_binary_procedures)[where], " procedures")
}
setkey(dt_binary_procedures, "eid")

dt_binary_classified <- merge(dt_binary_classified, dt_binary_procedures)
fwrite(dt_binary_classified, file="/well/lindgren/UKBIOBANK/dpalmer/superpopulation_assignments/BRaVa_phenotypes_with_superpopulation_labels.tsv", sep="\t")

for (i in 1:length(ICD_and_procedures)) {
	dt_binary_classified[[paste(ICD_and_procedures[i], "total")]] <- dt_binary_classified[[ICD_and_procedures[i]]]
	dt_binary_classified[[paste(ICD_and_procedures[i], "total")]][which(dt_binary_classified[[paste(ICD_and_procedures[i], "procedures")]] == 1)] <- 1
}

# Split by 1000G label and count
dt_counts <- dt_binary_classified %>% group_by(classification_strict) %>% summarise(across(c(union(dt_query_ICD$phenotype, dt_query_OPCS4$phenotype), paste(ICD_and_procedures, c("procedures", "total"))), sum, na.rm=TRUE))
dt_counts_t <- data.table::transpose(dt_counts, keep.names="phenotype", make.names="classification_strict")

# Finally, comma separate and combine the ICD9 and ICD10 codes together for inclusion on ICD_phecode tab of BRaVa_Nominate_Phenotypes spreadsheet.
fwrite(dt_counts_t, file="data/output/UKBB_case_counts.tsv", sep="\t")

dt_counts <- dt_binary_classified %>% summarise(across(c(union(dt_query_ICD$phenotype, dt_query_OPCS4$phenotype), paste(ICD_and_procedures, c("procedures", "total"))), sum, na.rm=TRUE))
dt_counts_t <- data.table::transpose(dt_counts, keep.names="phenotype")
names(dt_counts_t)[2] <- "count"
fwrite(dt_counts_t, file="data/output/UKBB_case_counts_total.tsv", sep="\t")

