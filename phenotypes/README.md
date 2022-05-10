# Determining case control status for ICD codes

In this short script, we include a set of functions to help munge data from the BRaVa_Nominate_Phenotypes google sheet, and use this to define regular expressions to then extract from a collection of phenotype data available for a set of samples.

`munge_BRaVa_ICD_proposals` downloads and extracts the ICD9/10 codes from the BRaVa_Nominate_Phenotypes google sheet ICD_Phecode tab. We clean it up, removing extra spaces and include a 'space pipe space' (` | `) delimiter to define the regular expression for case inclusions and control exclusions. 

`extract_ID_and_ICD_UKB` munges the relevent (ICD) columns of the UK Biobank data, and re-organises things so that each sample has a space delimited string of ICD10 and ICD9 codes associated to them, we then remove all columns except for `eid` (the sample ID), `ICD9_string` (a space delimited string of ICD9 codes) and `ICD10_string` (a space delimited string of ICD10 codes).

`extract_case_status` loops over all the ICD case inclusion and control exclusion definitions for each of the proposed ICD phenotypes and assigns each individual case (1), control (0), or missing (control exclusion) status.

Finally, we can use the output from `05_estimate_superpopulation.r` to assign superpopulation labels, by merging, and then use `dplyr` to extract case and control counts split by superpopulation.
