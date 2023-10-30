#!/bin/sh

# Load R module
module load R

# phenotypes=("Age_related_macular_degeneration"
#    "Alanine_transaminase"
#    "Alcohol_consumption_drinks_per_week"
#    "Aspartate_aminotransferase"
#    "Asthma"
#    "Atrial_Fibrillation"
#    "Benign_and_in_situ_intestinal_neoplasms"
#    "BMI"
#    "C_reactive_protein"
#    "Chronic_obstructive_pulmonary_disease"
#    "Chronic_Renal_Failure"
#    "Colon_and_rectum_cancer"
#    "Coronary_artery_disease"
#    "Gout"
#    "HDLC"
#    "Heart_Failure"
#    "Height"
#    "Hip_replacement"
#    "Hypertension"
#    "Inflammatory_bowel_disease"
#    "Inguinal_femoral_and_abdominal_hernia"
#    "Interstitial_lung_disease_and_pulmonary_sarcoidosis"
#    "LDLC"
#    "Non_rheumatic_valvular_heart_disease"
#    "Pancreatitis"
#    "Peptic_ulcer_disease"
#    "Peripheral_artery_disease"
#    "Psoriasis"
#    "Rheumatic_heart_disease"
#    "Rheumatoid_arthritis"
#    "Stroke"
#    "Total_cholesterol"
#    "Triglycerides"
#    "Type_2_diabetes"
#    "Urolithiasis"
#    "Varicose_Veins"
#    "Venous_Thromboembolism"
#    "WHR_adjusted_for_BMI")

phenotypes=("Age_related_macular_degeneration")
pops=("AFR" "AMR" "EAS" "EUR" "SAS")

for phenotype in "${phenotypes[@]}"
do
   for pop in "${pops[@]}"
   do
      # Submit sbatch, passing name as a parameter
      sbatch \
      --job-name=$phenotype \
      -o="$phenotype.out" \
      Rscript run_qq_and_manhattan.r \
      --phenotype=$phenotype \
      --population=$pop \
      --sex "both_sexes"
   done
done

# phenotypes=("Benign_and_in_situ_cervical_and_uterine_neoplasms"
#    "Breast_cancer"
#    "Cervical_cancer"
#    "Excess_frequent_and_irregular_menstrual_bleeding"
#    "Female_infertility"
#    "Maternal_hemorrhage")

# for phenotype in "${phenotypes[@]}"
# do
#    for pop in "${pops[@]}"
#    do
#       # Submit sbatch, passing name as a parameter
#       sbatch \
#       --job-name=$phenotype \
#       --o="$phenotype.out" \
#       Rscript run_qq_and_manhattan.r \
#       --phenotype=$phenotype \
#       --population=$pop \
#       --sex "XX"
#    done
# done









