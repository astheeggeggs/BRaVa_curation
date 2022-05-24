# install.packages("googlesheets4")
library(data.table)
library(dplyr)
library(ggplot2)
library(googlesheets4)

munge_BRaVa_ICD_proposals <- function() {
	# Download and extract the case and control ICD9/10 codes from the BRaVa nominate phenotypes file
	dt <- read_sheet("https://docs.google.com/spreadsheets/d/1YqdSyxf2OyoIYvLnDVj7NmbpebtppsgyJSq18gkVAWI/edit#gid=1716081249", sheet=3)
	cols <- c(
		"Description",
		"ICD_control_exclude (subjects with these ICD codes are excluded from controls )",
		"ICD_case_include (ICD codes to define cases)")
	dt <- dt[, cols, with=FALSE]
	names(dt) <- c("phenotype", "control_exclude", "case_include")
	dt <- data.table(dt)

	# Remove any extra spaces
	dt[, control_exclude := gsub("[[:space:]]", "", control_exclude)]
	dt[, case_include := gsub("[[:space:]]", "", case_include)]
	
	# Replace '.'s with '' (it's an equivalent encoding for ICD codes).
	dt[, control_exclude := gsub("\\.", "", control_exclude)]
	dt[, case_include := gsub("\\.", "", case_include)]

	# Replace commas for pipes to get the regular expression ready.
	dt[, control_exclude := gsub(",", " | ", control_exclude)]
	dt[, case_include := gsub(",", " | ", case_include)]

	# Ensure that the first and final character of the regular expression is a space.
	dt[, control_exclude := gsub("^(.*)$", " \\1 ", control_exclude)]
	dt[, case_include := gsub("^(.*)$", " \\1 ", case_include)]

	dt$control_exclude[which(is.na(dt$control_exclude))] <- ""
	dt$case_include[which(is.na(dt$case_include))] <- ""

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

extract_case_status <- function(dt_data, dt_query) {
	cols <- dt_query$phenotype
	dt_data[, (cols) := 0]
	for (i in 1:length(cols)) {
		# Add a new column to dt_data
		phenotype <- cols[i]
		cat(paste0(phenotype, "..."))
		if (dt_query$case_include[i] != "") {
			ICD_where <- grep(dt_query$case_include[i], dt_data$ICD10_string)
			cat(paste0("cases:", length(ICD_where), "..."))
			dt_data[[phenotype]][ICD_where] <- 1
			ICD_where <- grep(dt_query$case_include[i], dt_data$ICD9_string)
			cat(paste(length(ICD_where), "...\n"))
			dt_data[[phenotype]][ICD_where] <- 1
		}

		if (dt_query$control_exclude[i] != "") {
			ICD_where <- grep(dt_query$control_exclude[i], dt_data$ICD10_string)
			cat(paste0("control exclusions:", length(ICD_where), "..."))
			dt_data[[phenotype]][ICD_where] <- NA
			ICD_where <- grep(dt_query$control_exclude[i], dt_data$ICD9_string)
			cat(paste0(length(ICD_where), "...\n"))
			dt_data[[phenotype]][ICD_where] <- NA
		}
		cat(paste0("count: ", sum(dt_data[[phenotype]], na.rm=TRUE)))
	}
	return(dt_data)
}

dt_query <- munge_BRaVa_ICD_proposals()
fwrite(dt_query, file='data/BRaVa_ICD_proposals.tsv', sep='\t', quote=TRUE)
dt_query <- fread('data/BRaVa_ICD_proposals.tsv')
dt_data <- extract_ID_and_ICD_UKB()
dt_binary <- extract_case_status(dt_data, dt_query)
setkey(dt_binary, "eid")

# Merge with 1000G labels - this file is created using 05_estimate_superpopulation.r in the QC folder.
dt_classify <- fread("/well/lindgren/UKBIOBANK/dpalmer/superpopulation_assignments/superpopulation_labels.tsv")
dt_classify[, eid:=sample.ID]
dt_classify[, sample.ID:=NULL]
setkey(dt_classify, "eid")

dt_binary_classified <- merge(dt_binary, dt_classify)

# Split by 1000G label and count
dt_counts <- dt_binary_classified %>% group_by(classification_strict) %>% summarise(across(dt_query$phenotype, sum, na.rm=TRUE))

dt_counts_t <- data.table::transpose(dt_counts, keep.names="phenotype", make.names="classification_strict")

# Finally, comma separate and combine the ICD9 and ICD10 codes together for inclusion on ICD_phecode tab of BRaVa_Nominate_Phenotypes spreadsheet.
fwrite(dt_counts_t, file="data/output/UKBB_case_counts.tsv", sep="\t")

dt_counts <- dt_binary_classified %>% summarise(across(dt_query$phenotype, sum, na.rm=TRUE))
dt_counts_t <- data.table::transpose(dt_counts, keep.names="phenotype")
names(dt_counts_t)[2] <- "count"
fwrite(dt_counts_t, file="data/output/UKBB_case_counts_total.tsv", sep="\t")


