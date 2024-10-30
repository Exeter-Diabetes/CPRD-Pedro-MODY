
# Extract MODY cohort with new variable names for Pedro to calculate MODY scores


############################################################################################
  
# Setup
library(tidyverse)
library(aurum)
rm(list=ls())

cprd = CPRDData$new(cprdEnv = "diabetes-jun2024",cprdConf = "~/.aurum.yaml")

analysis = cprd$analysis("dpctn_final")


############################################################################################

# Get cohort info

cohort <- cohort %>% analysis$cached("cohort_24")


############################################################################################

# Define MODY cohort: patients diagnosed with a current Type 1 or Type 2 diagnosis, diagnosed aged 1-35
## Don't exclude if missing BMI/HbA1c before diagnosis (NB: missing BMI if no BMI in adulthood too) - different to main MODY cohort

pedro_mody_cohort <- cohort %>%
  filter(dm_diag_age>=1 & dm_diag_age<36 & (diabetes_type=="type 1" | diabetes_type=="type 2" | diabetes_type=="mixed; type 1" | diabetes_type=="mixed; type 2") & !is.na(diagnosis_date)) %>%
  mutate(c_peptide_ever=ifelse(!is.na(earliest_c_pep_ins_deficient) | !is.na(earliest_c_pep_ins_intermediate) | !is.na(earliest_c_pep_ins_normal), 1L, 0L),
         abs_ever=ifelse(!is.na(earliest_negative_gad) | !is.na(earliest_positive_gad) | !is.na(earliest_negative_ia2) | !is.na(earliest_positive_ia2), 1L, 0L),
         hba1c_mmolmol=ifelse(hba1cdate>=diagnosis_date, hba1c, NA),
         hba1c=(0.09148*hba1c_mmolmol)+2.152,
         age_at_bmi=datediff(bmidate, dob)/365.25,
         bmi=ifelse(bmidate>=diagnosis_date & age_at_bmi>=18, bmi, NA),
         sex=gender-1,
         insoroha=ifelse(current_oha==1 | current_insulin==1, 1L, 0L),
         
         which_equation=ifelse((current_insulin==1 & current_mfn==0 & current_su==0 & current_dpp4==0 & current_sglt2==0 & current_gipglp1==0 & current_glp1==0 & current_ldSema==0 & current_hdSema==0 & current_oSema==0 & current_tzd==0) | ((diabetes_type=="type 1" | diabetes_type=="mixed; type 1") & dm_diag_age<10 & hba1c_mmolmol>58) | ((diabetes_type=="type 1" | diabetes_type=="mixed; type 1") & dm_diag_age<25 & current_insulin==1 & current_mfn==1 & current_su==0 & current_dpp4==0 & current_sglt2==0 & current_gipglp1==0 & current_glp1==0 & current_ldSema==0 & current_hdSema==0 & current_oSema==0 & current_tzd==0), "t1", "t2")) %>%
  
  
  select(patid, ethnicity_5cat, ethnicity_16cat, c_peptide_ever, abs_ever, agedx=dm_diag_age, sex, hba1c, bmi, pardm=fh_diabetes, insoroha, agerec=age_at_index, diabetes_type, which_equation) %>%
  
  analysis$cached("pedro_mody_cohort_24", unique_indexes="patid")

pedro_mody_cohort %>% count()
#106230


#c-peptide and ab counts
pedro_mody_cohort %>% filter(c_peptide_ever==1 & abs_ever==1) %>% count()
#387

387/106230
#0.4%

