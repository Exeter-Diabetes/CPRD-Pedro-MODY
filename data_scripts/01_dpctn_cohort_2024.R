
# Create new table with all IDs in download who are registered on 01/02/2020 and adults at this date

# Define diagnosis date based on earliest diabetes code: exclude those diagnosed between -30 and +90 days of registration or >50 years of age

# Define diabetes type based on (latest) type codes - and add in variables relating to this

# Pull in variables for MODY and T1D/T2D calculator: current BMI, HbA1c, total cholesterol, HDL, triglycerides, current treatment (and whether have ins/OHA script or DPP4i/GLP1/SU/TZD ever), family history of diabetes, earliest insulin and OHA

# Pull in other variables of interest: GAD and IA2 antibodies (ever prior to index date), C-peptide (ever prior to index date), and whether are non-English speaking / have English as a second language

# NB: don't have linkage data at this stage: use cprd_ddate for death, will have higher missing ethnicity, and won't have IMD info


# All is using latest codelists (OHA table still in old format but doesn't affect results)


############################################################################################

# Setup
library(tidyverse)
library(aurum)
library(EHRBiomarkr)
rm(list=ls())

cprd = CPRDData$new(cprdEnv = "diabetes-jun2024",cprdConf = "~/.aurum.yaml")
codesets = cprd$codesets()
codes = codesets$getAllCodeSetVersion(v = "01/06/2024")

analysis = cprd$analysis("dpctn_final")


############################################################################################

# Set index date

index_date <- as.Date("2024-01-01")


############################################################################################

# Find IDs of those registered on 01/02/2020

# Include general patient variables as per all_diabetes_cohort:
## gender
## DOB (derived as per: https://github.com/Exeter-Diabetes/CPRD-Codelists/blob/main/readme.md#general-notes-on-implementation, cached for all IDs in 'diagnosis_date_dob' table from all_patid_dob script)
## pracid
## prac_region
## ethnicity 5-category and 16-category (derived as per: https://github.com/Exeter-Diabetes/CPRD-Codelists#ethnicity, from GP data only)
## regstartdate	
## gp_record_end (earliest of last collection date from practice, deregistration and 30/06/2024 (latest date in records))	
## death_date	('cprddeathdate')
## Calculate age at index date


# Just want those registered on 01/02/2020 (the index date)

analysis = cprd$analysis("diabetes_cohort")
dob <- dob %>% analysis$cached("dob")

analysis = cprd$analysis("all_patid")
ethnicity <- ethnicity %>% analysis$cached("ethnicity")


analysis = cprd$analysis("dpctn_final")

cohort <- cprd$tables$patient %>%
  select(patid, gender, pracid, regstartdate, regenddate, cprd_ddate) %>%
  left_join((dob %>% select(patid, dob)), by="patid") %>%
  left_join((cprd$tables$practice %>% select(pracid, prac_region=region, lcd)), by="pracid") %>%
  left_join((ethnicity %>% select(patid, ethnicity_5cat, ethnicity_16cat)), by="patid") %>%
  
  mutate(gp_record_end=pmin(if_else(is.na(lcd), as.Date("2024-06-30"), lcd),
                            if_else(is.na(regenddate), as.Date("2024-06-30"), regenddate),
                            as.Date("2024-06-30"), na.rm=TRUE),
         
         death_date=cprd_ddate,
         
         age_at_index=round(datediff(index_date, dob)/365.25, 1)) %>%
  
  select(patid, gender, dob, age_at_index, pracid, prac_region, ethnicity_5cat, ethnicity_16cat, regstartdate, gp_record_end, death_date) %>%
  
  filter(regstartdate<=index_date & !(!is.na(death_date) & death_date<index_date) & !(!is.na(gp_record_end) & gp_record_end<index_date)) %>%
  
  analysis$cached("cohort_24_interim_1", unique_indexes="patid")

cohort %>% count()
# 1,323,727


# Remove patients from practices which may have duplicated data as per CPRD's guidance

analysis = cprd$analysis("diabetes_cohort")

practice_exclusion_ids <- cprd$tables$patient %>% 
  filter(pracid == "20024" | pracid == "20036" |pracid == "20091" |pracid == "20171" | pracid == "20178" |pracid == "20202" | pracid == "20254" | pracid == "20389" |pracid == "20430" |pracid == "20452" |
           pracid == "20469" | pracid == "20487" | pracid == "20552" | pracid == "20554" | pracid == "20640" | pracid == "20717" | pracid == "20734" | pracid == "20737" | pracid == "20740" | pracid == "20790" |
           pracid == "20803" | pracid == "20822" | pracid == "20868" | pracid == "20912" | pracid == "20996" | pracid == "21001" | pracid == "21015" | pracid == "21078" | pracid == "21112" | pracid == "21118" |
           pracid == "21172" | pracid == "21173" | pracid == "21277" | pracid == "21281" | pracid == "21331" | pracid == "21334" | pracid == "21390" | pracid == "21430" | pracid == "21444" | pracid == "21451" |
           pracid == "21529" | pracid == "21553" | pracid == "21558" | pracid == "21585") %>%
  distinct(patid) %>%
  analysis$cached("practice_exclusion_ids")


analysis = cprd$analysis("dpctn_final")

cohort <- cohort %>%
  anti_join(practice_exclusion_ids, by="patid") %>%
  analysis$cached("cohort_24_interim_2", unique_indexes="patid")

cohort %>% count()
# 1,323,727



# Just keep men and women

cohort <- cohort %>%
  filter(gender==1 | gender==2) %>%
  analysis$cached("cohort_24_interim_3", unique_indexes="patid")

cohort %>% count()
# 1,323,701



# Just keep those which are adults

cohort <- cohort %>%
  filter(age_at_index>=18) %>%
  analysis$cached("cohort_24_interim_4", unique_indexes="patid")

cohort %>% count()
# 1,304,837


############################################################################################

# Define diabetes type (NB: not removing codes before DOB for this)
# Have changed script in Aug 2024 so only use dpctn IDs as quicker - and have renamed tables accordingly so not 'all_patid'

## Find all diabetes codes prior to/at index date

analysis = cprd$analysis("all_patid")

raw_diabetes_medcodes <- cprd$tables$observation %>%
  inner_join(codes$all_diabetes, by="medcodeid") %>%
  analysis$cached("raw_diabetes_medcodes", indexes=c("patid", "obsdate", "all_diabetes_cat"))


analysis = cprd$analysis("dpctn_final")

dm_codes <- cohort %>%
  left_join(raw_diabetes_medcodes, by="patid") %>%
  select(patid, date=obsdate, enterdate, category=all_diabetes_cat) %>%
  filter(date<=index_date) %>%
  analysis$cached("dm_codes_24", indexes=c("patid", "date", "category"))

## Check how have diabetes codes before index date in download - won't be everyone
dm_codes %>% distinct(patid) %>% count()
#1,278,300 - have lost 26,537


## Find code counts for each diabetes type

dm_code_counts <- dm_codes %>%
  group_by(patid, category) %>%
  summarise(count=n()) %>%
  ungroup() %>%
  pivot_wider(id_cols=patid, names_from=category, values_from=count, values_fill=list(count=0)) %>%
  analysis$cached("dm_code_counts_24", indexes="patid")

## Are insulin receptor abs now
## Now have ambiguous pregnancy codes - fine, don't use for classification as not clear which category - will be unclassified if only have these codes

## Find latest type code, excluding unspecified and gestational, and keep date
 
latest_type_code <- dm_codes %>%
  filter(category!="unspecified" & category!="gestational" & category!="gestational history" & category!="pregnancy") %>%
  group_by(patid) %>%
  mutate(most_recent_date=max(date, na.rm=TRUE),
         days_since_type_code=datediff(index_date, most_recent_date)) %>%
  filter(date==most_recent_date) %>%
  ungroup() %>%
  group_by(patid, days_since_type_code) %>%
  summarise(provisional_diabetes_type=sql("group_concat(distinct category order by category separator ' & ')")) %>%
  ungroup() %>%
  analysis$cached("latest_type_code_24", indexes="patid")


## Determine diabetes type and add to cohort table

diabetes_type <- dm_code_counts %>%
  left_join(latest_type_code, by="patid") %>%
  mutate(diabetes_type=case_when(
    `type 1`==0 & `type 2`==0 & gestational==0 & `gestational history`==0 & malnutrition==0 & mody==0 & `other type not specified`==0 & `other/unspec genetic inc syndromic`==0 & secondary==0 & "insulin receptors abs"==0 ~ "unspecified",
    
    `type 1`>0 & `type 2`==0 & gestational==0 & `gestational history`==0 & malnutrition==0 & mody==0 & `other type not specified`==0 & `other/unspec genetic inc syndromic`==0 & secondary==0 & "insulin receptors abs"==0 ~ "type 1",
    
    `type 1`==0 & `type 2`>0 & gestational==0 & `gestational history`==0 & malnutrition==0 & mody==0 & `other type not specified`==0 & `other/unspec genetic inc syndromic`==0 & secondary==0 & "insulin receptors abs"==0 ~ "type 2",
    
    `type 1`==0 & `type 2`==0 & (gestational>0 | `gestational history`>0) & malnutrition==0 & mody==0 & `other type not specified`==0 & `other/unspec genetic inc syndromic`==0 & secondary==0 & "insulin receptors abs"==0 ~ "gestational",
    
    `type 1`==0 & `type 2`==0 & gestational==0 & `gestational history`==0 & malnutrition==0 & mody==0 & `other type not specified`==0 & `other/unspec genetic inc syndromic`==0 & secondary==0 & "insulin receptors abs">0 ~ "insulin receptor abs",
    
    `type 1`==0 & `type 2`==0 & gestational==0 & `gestational history`==0 & malnutrition>0 & mody==0 & `other type not specified`==0 & `other/unspec genetic inc syndromic`==0 & secondary==0 & "insulin receptors abs"==0 ~ "malnutrition",
    
    `type 1`==0 & `type 2`==0 & gestational==0 & `gestational history`==0 & malnutrition==0 & mody>0 & `other type not specified`==0 & `other/unspec genetic inc syndromic`==0 & secondary==0 & "insulin receptors abs"==0 ~ "mody",
    
    `type 1`==0 & `type 2`==0 & gestational==0 & `gestational history`==0 & malnutrition==0 & mody==0 & `other type not specified`>0 & `other/unspec genetic inc syndromic`>0 & secondary==0 & "insulin receptors abs"==0 ~ "other type not specified",
    
    `type 1`==0 & `type 2`==0 & gestational==0 & `gestational history`==0 & malnutrition==0 & mody==0 & `other type not specified`==0 & `other/unspec genetic inc syndromic`>0 & secondary==0 & "insulin receptors abs"==0 ~ "other/unspec genetic inc syndromic",
    
    `type 1`==0 & `type 2`==0 & gestational==0 & `gestational history`==0 & malnutrition==0 & mody==0 & `other type not specified`==0 & `other/unspec genetic inc syndromic`==0 & secondary>0 & "insulin receptors abs"==0 ~ "secondary"),
    
    diabetes_type=ifelse(!is.na(diabetes_type), diabetes_type, paste("mixed;", provisional_diabetes_type)),
    
    type1_code_count= `type 1`,
    type2_code_count= `type 2`,
    mody_code_count=mody
    
    ) %>%
  
  select(patid, diabetes_type, type1_code_count, type2_code_count, mody_code_count, days_since_type_code) %>%
  
  analysis$cached("diabetes_type_24", unique_indexes="patid")
    

cohort <- cohort %>%
  inner_join(diabetes_type, by="patid") %>%
  analysis$cached("cohort_24_interim_5", unique_indexes="patid")

cohort %>% count()
#1,278,300


############################################################################################

# Define diagnosis dates
## Remove if before DOB
## Exclude codes if Type 2 and in year in birth

## Check how have diabetes codes before DOB only, or with T2 and codes in YOB only
dm_codes %>% inner_join(cprd$tables$validDateLookup, by="patid") %>% filter(date>=min_dob) %>% distinct(patid) %>% count()
#1,278,300 - lose noone with codes before DOB only

diagnosis_dates <- dm_codes %>%
  inner_join(cprd$tables$validDateLookup, by="patid") %>%
  filter(date>=min_dob) %>%
  inner_join(cohort, by="patid") %>%
  filter(!(diabetes_type=="type 2" & year(date)==year(dob))) %>%
  group_by(patid) %>%
  summarise(diagnosis_date=min(date, na.rm=TRUE)) %>%
  ungroup() %>%
  analysis$cached("diagnosis_dates_24", unique_indexes="patid")

diagnosis_dates %>% count()
#still 1,278,300, so also no-one with T2 and codes in YOB only


# Add to cohort table (no restrictions on diagnosis age)
cohort <- cohort %>%
  inner_join(diagnosis_dates, by="patid") %>%
  mutate(dm_diag_age=round((datediff(diagnosis_date, dob))/365.25, 1)) %>%
  analysis$cached("cohort_24_interim_6", unique_indexes="patid")

cohort %>% count()
#1,278,300


# Set diagnosis date to missing where within -30 to +90 days of registration start
## Keep age at diagnosis

cohort <- cohort %>%
  mutate(diagnosis_date=if_else(datediff(diagnosis_date, regstartdate)< -30 | datediff(diagnosis_date, regstartdate)>90, diagnosis_date, as.Date(NA))) %>%
  analysis$cached("cohort_24_interim_7", unique_indexes="patid")


############################################################################################

# Add in biomarkers

## Don't set any limit on how far back to go - will remove those before diagnosis

# Clean biomarkers:
## Only keep those within acceptable value limits
## Only keep those with valid unit codes (numunitid)
## If multiple values on the same day, take mean
## Remove those with invalid dates (before DOB or after LCD/death/deregistration)

# Find baseline values
## Use closest date to index date as long as prior to this


biomarkers <- c("bmi", "hdl", "triglyceride", "totalcholesterol", "hba1c")


# Pull out all raw biomarker values and cache

analysis = cprd$analysis("all_patid")

for (i in biomarkers) {
  
  print(i)
  
  raw_tablename <- paste0("raw_", i, "_medcodes")
  
  data <- cprd$tables$observation %>%
    inner_join(codes[[i]], by="medcodeid") %>%
    analysis$cached(raw_tablename, indexes=c("patid", "obsdate", "testvalue", "numunitid"))
  
  assign(raw_tablename, data)
  
}


# Clean biomarkers:
## Only keep those within acceptable value limits
## Only keep those with valid unit codes (numunitid)
## If multiple values on the same day, take mean
## Remove those with invalid dates (before min DOB or after LCD/death/deregistration)

## HbA1c only: remove before 1990, and convert all values to mmol/mol

analysis = cprd$analysis("all_patid")

for (i in biomarkers) {
  
  print(i)
  
  raw_tablename <- paste0("raw_", i, "_medcodes")
  clean_tablename <- paste0("clean_", i, "_medcodes")
  
  data <- get(raw_tablename)
  
  if (i=="hba1c") {
    
    data <- data %>%
      filter(year(obsdate)>=1990) %>%
      mutate(testvalue=ifelse(testvalue<=20, ((testvalue-2.152)/0.09148), testvalue))
        
  }
    
  data <- data %>%
    
    clean_biomarker_values(testvalue, i) %>%
    clean_biomarker_units(numunitid, i) %>%
    
    group_by(patid,obsdate) %>%
    summarise(testvalue=mean(testvalue, na.rm=TRUE)) %>%
    ungroup() %>%
    
    inner_join(cprd$tables$validDateLookup, by="patid") %>%
    filter(obsdate>=min_dob & obsdate<=gp_ons_end_date) %>%
    
    select(patid, date=obsdate, testvalue) %>%
    
    analysis$cached(clean_tablename, indexes=c("patid", "date", "testvalue"))
  
  assign(clean_tablename, data)
  
}


# For each biomarker, find baseline value at index date
## Use closest date to index date as long as prior to this - weight and height don't need to be on same day as height doesn't change

analysis = cprd$analysis("dpctn_final")


for (i in biomarkers) {
  
  print(i)
  
  clean_tablename <- paste0("clean_", i, "_medcodes")
  biomarker_date_variable <- paste0(i, "date")
  biomarker_indexdiff_variable <- paste0(i, "indexdiff")
  
  data <- get(clean_tablename) %>%
    mutate(indexdatediff=datediff(date, index_date)) %>%
    filter(indexdatediff<=0) %>%
    group_by(patid) %>%
    mutate(min_timediff=min(abs(indexdatediff), na.rm=TRUE)) %>%
    filter(abs(indexdatediff)==min_timediff) %>%
    ungroup() %>%
    select(patid, testvalue, date, indexdatediff) %>%
    rename({{i}}:=testvalue,
           {{biomarker_date_variable}}:=date,
           {{biomarker_indexdiff_variable}}:=indexdatediff)
  
  cohort <- cohort %>%
    left_join(data, by="patid")
  
}

cohort <- cohort %>% analysis$cached("cohort_24_interim_8", unique_indexes="patid")


############################################################################################

# Add in GAD and IA2 antibodies

analysis = cprd$analysis("all_patid")

# GAD
raw_gad_medcodes <- cprd$tables$observation %>%
  inner_join(codes$gad, by="medcodeid") %>%
  analysis$cached("raw_gad_medcodes", indexes=c("patid", "obsdate", "testvalue", "numunitid"))
raw_gad_medcodes %>% count()
#11416

# IA2
raw_ia2_medcodes <- cprd$tables$observation %>%
  inner_join(codes$ia2, by="medcodeid") %>%
  analysis$cached("raw_ia2_medcodes", indexes=c("patid", "obsdate", "testvalue", "numunitid"))
raw_ia2_medcodes %>% count()
#9042


print(raw_gad_medcodes %>% group_by(numunitid) %>% summarise(count=n()), n=100)
print(raw_ia2_medcodes %>% group_by(numunitid) %>% summarise(count=n()), n=100)

# Assume those with missing units are in U/mL
## GAD: majority are missing (5972), then rest are (i)U/mL / k(i)U/L (equivalent) / 'units'; 1 mmol/L, 1 mu/L, 1 u/L, and 2 n/ml - exclude
## IA2: vast majority are missing (8317), then some (i)U/mL / k(i)U/L (equivalent) / 'units'; 3 'unknown UoM', 2 u/L, 1 'titre' - exclude
## Use >=11 U/mL as positive for GAD; >=7.5 U/mL for IA2 as per: https://www.exeterlaboratory.com/

# Assume 0 values are 0, not missing

analysis = cprd$analysis("dpctn_final")

### GAD

gad <- raw_gad_medcodes %>%
  filter(obsdate<=index_date & (is.na(numunitid) | (numunitid!=229 & numunitid!=11589 & numunitid!=218 & numunitid!=276))) %>%
  mutate(result=ifelse(!is.na(testvalue) & testvalue<11, "negative",
                    ifelse(!is.na(testvalue) & testvalue>=11, "positive", NA))) %>%
  filter(!is.na(result)) %>%
  distinct(patid, obsdate, result)

earliest_gad <- gad %>%
  group_by(patid, result) %>%
  mutate(earliest_gad=min(obsdate, na.rm=TRUE)) %>%
  filter(obsdate==earliest_gad) %>%
  ungroup() %>%
  pivot_wider(id_cols="patid", names_from=result, values_from=obsdate)

latest_gad <- gad %>%
  group_by(patid, result) %>%
  mutate(latest_gad=max(obsdate, na.rm=TRUE)) %>%
  filter(obsdate==latest_gad) %>%
  ungroup() %>%
  pivot_wider(id_cols="patid", names_from=result, values_from=obsdate)

gad <- earliest_gad %>%
  rename(earliest_negative_gad=negative, earliest_positive_gad=positive) %>%
  left_join((latest_gad %>% rename(latest_negative_gad=negative, latest_positive_gad=positive)), by="patid") %>%
  select(patid, earliest_negative_gad, latest_negative_gad, earliest_positive_gad, latest_positive_gad) %>%
  analysis$cached("gad_24", unique_indexes="patid")


### IA2

ia2 <- raw_ia2_medcodes %>%
  filter(obsdate<=index_date & (is.na(numunitid) | (numunitid!=276 & numunitid!=292 & numunitid!=274))) %>%
  mutate(result=ifelse(!is.na(testvalue) & testvalue<7.5, "negative",
                       ifelse(!is.na(testvalue) & testvalue>=7.5, "positive", NA))) %>%
  filter(!is.na(result)) %>%
  distinct(patid, obsdate, result)

earliest_ia2 <- ia2 %>%
  group_by(patid, result) %>%
  mutate(earliest_ia2=min(obsdate, na.rm=TRUE)) %>%
  filter(obsdate==earliest_ia2) %>%
  ungroup() %>%
  pivot_wider(id_cols="patid", names_from=result, values_from=obsdate)

latest_ia2 <- ia2 %>%
  group_by(patid, result) %>%
  mutate(latest_ia2=max(obsdate, na.rm=TRUE)) %>%
  filter(obsdate==latest_ia2) %>%
  ungroup() %>%
  pivot_wider(id_cols="patid", names_from=result, values_from=obsdate)

ia2 <- earliest_ia2 %>%
  rename(earliest_negative_ia2=negative, earliest_positive_ia2=positive) %>%
  left_join((latest_ia2 %>% rename(latest_negative_ia2=negative, latest_positive_ia2=positive)), by="patid") %>%
  select(patid, earliest_negative_ia2, latest_negative_ia2, earliest_positive_ia2, latest_positive_ia2) %>%
  analysis$cached("ia2_24", unique_indexes="patid")


cohort <- cohort %>%
  left_join(gad, by="patid") %>%
  left_join(ia2, by="patid") %>%
  analysis$cached("cohort_24_interim_9", unique_indexes="patid")


############################################################################################

# Add in C-peptide

analysis = cprd$analysis("all_patid")

raw_c_peptide_medcodes <- cprd$tables$observation %>%
  inner_join(codes$c_peptide, by="medcodeid") %>%
  analysis$cached("raw_c_peptide_medcodes", indexes=c("patid", "obsdate", "testvalue", "numunitid"))
raw_c_peptide_medcodes %>% count()
#15,050


# Look at units

raw_c_peptide_medcodes %>% filter(c_peptide_cat=="ucpcr") %>% group_by(numunitid) %>% summarise(count=n())
## 1154 are 2434 (nM/mM crea), 891 are 899 (nmol/mmol), 185 are missing, 43 are 959 (umol/mol), 1 is 1744 (nmol/mmol(creat)) - all good and all the same unit
## 6 are 235 (nmol/L), 1 each of 218 (mmol/L) and 233 (ng/mL) - exclude

raw_c_peptide_medcodes %>% filter(c_peptide_cat=="blood") %>% group_by(numunitid) %>% summarise(count=n())
## 7909 are 256 (pmol/L), 3923 are missing, 577 are 235 (nmol/L), 2 are 218 (mmol/L), 2 are 9244 (pmol/ml), 1 is 2339 (pm/L) - all good but need to convert all to pmol/L
### Convert 235 to pmol/L (divide by 1000)
### Convert 218 to pmol/L (divide by 10^6)
### Convert 9244 to pmol/L (multiply by 1000)
## 293 are 283 (ug/L), 36 are 899 (nmol/mmol), 24 are 277 (U/mL), 1 is 229 (mu/L) - exclude


# Clean and define insulin status
## Use thresholds of <0.2 and >=0.6 (UCPCR) / <200 and >=600 (blood) to define insulin status (https://www.exeterlaboratory.com/test/c-peptide-urine/)
### No indication as to whether blood samples are fasted or not so assume not
## Assume 0 values are 0, not missing

analysis = cprd$analysis("dpctn_final")

processed_c_peptide <- raw_c_peptide_medcodes %>%
  filter(obsdate<=index_date & !is.na(testvalue) & !(c_peptide_cat=="ucpcr" & (numunitid==235 | numunitid==218 | numunitid==233)) & !(c_peptide_cat=="blood" & (numunitid==283 | numunitid==899 | numunitid==277 | numunitid==229))) %>%
  mutate(new_testvalue=ifelse(numunitid==235, testvalue/1000,
                              ifelse(numunitid==218, testvalue/1000000,
                                     ifelse(numunitid==9244, testvalue*1000, testvalue)))) %>%
  mutate(c_pep_insulin=ifelse(c_peptide_cat=="ucpcr" & new_testvalue<0.2, "c_pep_ins_deficient",
                              ifelse(c_peptide_cat=="ucpcr" & new_testvalue>=0.2 & testvalue<0.6, "c_pep_ins_intermediate",
                                     ifelse(c_peptide_cat=="ucpcr" & new_testvalue>=0.6, "c_pep_ins_normal",
                                            ifelse(c_peptide_cat=="blood" & new_testvalue<200, "c_pep_ins_deficient",
                                                   ifelse(c_peptide_cat=="blood" & new_testvalue>=200 & new_testvalue<600, "c_pep_ins_intermediate",
                                                          ifelse(c_peptide_cat=="blood" & new_testvalue>=600, "c_pep_ins_normal", NA))))))) %>%
  select(patid, date=obsdate, testvalue=new_testvalue, c_pep_cat=c_peptide_cat, c_pep_insulin) %>%
  distinct() %>%
  analysis$cached("processed_c_peptide_24", indexes=c("patid", "date", "c_pep_insulin"))



# check converted ones
test <- raw_c_peptide_medcodes %>%
  
  filter(obsdate<=index_date & !is.na(testvalue) & c_peptide_cat=="blood" & (numunitid==235 | numunitid==218 | numunitid==9244)) %>%
  
  mutate(original_units=ifelse(c_peptide_cat=="ucpcr", "urine_nmmm",
                               ifelse(numunitid==235, "nmol/L",
                                      ifelse(numunitid==218, "mmol/L",
                                             ifelse(numunitid==9244, "pmol/mL", "pmol/L"))))) %>%
  
  mutate(new_testvalue=ifelse(numunitid==235, testvalue/1000,
                              ifelse(numunitid==218, testvalue/1000000,
                                     ifelse(numunitid==9244, testvalue*1000, testvalue))),
         
         c_pep_insulin_original=ifelse(testvalue<200, "c_pep_ins_deficient",
                                       ifelse(testvalue>=200 & testvalue<600, "c_pep_ins_intermediate", "c_pep_ins_normal")),
         
         c_pep_insulin_new=ifelse(new_testvalue<200, "c_pep_ins_deficient",
                                  ifelse(testvalue>=200 & testvalue<600, "c_pep_ins_intermediate", "c_pep_ins_normal"))) %>%
  
  select(patid, date=obsdate, testvalue, c_pep_insulin_original, original_units, new_testvalue, c_pep_insulin_new) %>%
  collect()
#576 values

test %>% filter(c_pep_insulin_original!=c_pep_insulin_new) %>% count()
#5: converting units only affects 3 results anyway



# Add earliest and latest result for each category to table

earliest_c_pep <- processed_c_peptide %>%
  group_by(patid, c_pep_insulin) %>%
  mutate(earliest=min(date, na.rm=TRUE)) %>%
  filter(date==earliest) %>%
  ungroup() %>%
  pivot_wider(id_cols="patid", names_from=c_pep_insulin, values_from=date)

latest_c_pep <- processed_c_peptide %>%
  group_by(patid, c_pep_insulin) %>%
  mutate(latest=max(date, na.rm=TRUE)) %>%
  filter(date==latest) %>%
  ungroup() %>%
  pivot_wider(id_cols="patid", names_from=c_pep_insulin, values_from=date)

c_pep <- earliest_c_pep %>%
  rename(earliest_c_pep_ins_deficient=c_pep_ins_deficient, earliest_c_pep_ins_intermediate=c_pep_ins_intermediate, earliest_c_pep_ins_normal=c_pep_ins_normal) %>%
  left_join((latest_c_pep %>% rename(latest_c_pep_ins_deficient=c_pep_ins_deficient, latest_c_pep_ins_intermediate=c_pep_ins_intermediate, latest_c_pep_ins_normal=c_pep_ins_normal)), by="patid") %>%
  select(patid, earliest_c_pep_ins_deficient, latest_c_pep_ins_deficient, earliest_c_pep_ins_intermediate, latest_c_pep_ins_intermediate, earliest_c_pep_ins_normal, latest_c_pep_ins_normal) %>%
  analysis$cached("c_pep_24", unique_indexes="patid")


cohort <- cohort %>%
  left_join(c_pep, by="patid") %>%
  analysis$cached("cohort_24_interim_10", unique_indexes="patid")


############################################################################################

# Add in current treatment (last year)

# Get clean OHA and insulin scripts

analysis = cprd$analysis("all_patid")

## All OHA scripts
raw_oha_prodcodes <- cprd$tables$drugIssue %>%
  inner_join(cprd$tables$ohaLookup, by="prodcodeid") %>%
  analysis$cached("raw_oha_prodcodes", indexes=c("patid", "issuedate"))

clean_oha_prodcodes <- raw_oha_prodcodes %>%
  inner_join(cprd$tables$validDateLookup, by="patid") %>%
  filter(issuedate>=min_dob & issuedate<=gp_end_date) %>%
  select(patid, date=issuedate, dosageid, quantity, quantunitid, duration, Acarbose, DPP4, GIPGLP1, Glinide, GLP1, ldSema, hdSema, oSema, sema_query, MFN, SGLT2, SU, TZD, INS) %>%
  analysis$cached("clean_oha_prodcodes", indexes=c("patid", "date", "Acarbose", "DPP4", "GIPGLP1", "Glinide", "GLP1", "ldSema", "hdSema", "oSema", "sema_query", "MFN", "SGLT2", "SU", "TZD", "INS"))

## All insulin scripts
raw_insulin_prodcodes <- cprd$tables$drugIssue %>%
  inner_join(codes$insulin, by="prodcodeid") %>%
  analysis$cached("raw_insulin_prodcodes", indexes=c("patid", "issuedate"))

clean_insulin_prodcodes <- raw_insulin_prodcodes %>%
  inner_join(cprd$tables$validDateLookup, by="patid") %>%
  filter(issuedate>=min_dob & issuedate<=gp_end_date) %>%
  select(patid, date=issuedate, dosageid, quantity, quantunitid, duration, insulin_cat) %>%
  analysis$cached("clean_insulin_prodcodes", indexes=c("patid", "date", "insulin_cat"))

analysis = cprd$analysis("dpctn_final")


# Do drug classes separately and then combine (ignore Acarbose)

current_meds <- clean_oha_prodcodes %>%
  filter(date<=index_date & datediff(index_date, date)<=366) %>%
  group_by(patid) %>%
  summarise(current_dpp4=any(DPP4),
            current_gipglp1=any(GIPGLP1),
            current_glinide=any(Glinide),
            current_glp1=any(GLP1),
            current_ldSema=any(ldSema),
            current_hdSema=any(hdSema),
            current_oSema=any(oSema),
            current_mfn=any(MFN),
            current_sglt2=any(SGLT2),
            current_su=any(SU),
            current_tzd=any(TZD)) %>%
  ungroup() %>%
  analysis$cached("current_meds_24", unique_indexes="patid")


## Add insulin from insulin scripts (all prodcodes with OHA-insulin mixes are also in insulin)

current_insulin <- clean_insulin_prodcodes %>%
  filter(date<=index_date & datediff(index_date, date)<=366) %>%
  distinct(patid) %>%
  mutate(current_insulin=1L) %>%
  analysis$cached("current_insulin_24", unique_indexes="patid")


## Bolus/mix insulin

current_bolusmix_insulin <- clean_insulin_prodcodes %>%
  filter(date<=index_date & datediff(index_date, date)<=366 & (insulin_cat=="Bolus insulin" | insulin_cat=="Mix insulin")) %>%
  distinct(patid) %>%
  mutate(current_bolusmix_insulin=1L) %>%
  analysis$cached("current_bolusmix_insulin_24", unique_indexes="patid")



# Combine and make current_oha for any OHA

cohort <- cohort %>%
  left_join(current_meds, by="patid") %>%
  left_join(current_insulin, by="patid") %>%
  left_join(current_bolusmix_insulin, by="patid") %>%
  mutate(across(c("current_dpp4",
                  "current_glinide",
                  "current_gipglp1",
                  "current_glp1",
                  "current_ldSema",
                  "current_hdSema",
                  "current_oSema",
                  "current_mfn",
                  "current_sglt2",
                  "current_su",
                  "current_tzd",
                  "current_insulin",
                  "current_bolusmix_insulin"), coalesce, 0L)) %>%
  mutate(current_oha=ifelse(current_dpp4==1 | current_gipglp1==1 | current_glinide==1 | current_glp1==1 | current_ldSema==1 | current_hdSema==1 | current_oSema==1 | current_mfn==1 | current_sglt2==1 | current_su==1 | current_tzd==1, 1L, 0L)) %>%
  analysis$cached("cohort_24_interim_11", unique_indexes="patid")


############################################################################################

# Add in family history of diabetes (cleaned: includes 99% of raw occurrences so no difference)

# For people with positive and negative codes:
## If all negative codes are earlier than positive codes, fine - use positive
## Otherwise, treat as missing

analysis = cprd$analysis("all_patid")

## Raw FH of diabetes
raw_fh_diabetes_medcodes <- cprd$tables$observation %>%
  inner_join(codes$fh_diabetes, by="medcodeid") %>%
  analysis$cached("raw_fh_diabetes_medcodes", indexes=c("patid", "obsdate"))

## Clean FH of diabetes
clean_fh_diabetes_medcodes <- raw_fh_diabetes_medcodes %>%
  inner_join(cprd$tables$validDateLookup, by="patid") %>%
  filter(obsdate>=min_dob & obsdate<=gp_end_date) %>%
  select(patid, date=obsdate, fh_diabetes_cat) %>%
  analysis$cached("clean_fh_diabetes_medcodes", indexes=c("patid", "date", "fh_diabetes_cat"))


analysis = cprd$analysis("dpctn_final")

fh_code_types <- clean_fh_diabetes_medcodes %>%
  filter(fh_diabetes_cat!="positive - sibling" & fh_diabetes_cat!="positive - child" & fh_diabetes_cat!="positive - gestational" & date<=index_date) %>%
  mutate(fh_diabetes_cat=ifelse(fh_diabetes_cat=="negative", "negative", "positive")) %>%
  group_by(patid, fh_diabetes_cat) %>%
  summarise(earliest_date=min(date, na.rm=TRUE),
            latest_date=max(date, na.rm=TRUE)) %>%
  ungroup() %>%
  group_by(patid) %>%
  pivot_wider(id_cols=patid, names_from = c(fh_diabetes_cat), names_glue = "{fh_diabetes_cat}_{.value}", values_from=c(earliest_date, latest_date)) %>%
  ungroup() %>%
  analysis$cached("fh_code_types_24", unique_indexes="patid")

final_fh <- fh_code_types %>%
  mutate(fh_diabetes=ifelse(is.na(positive_earliest_date), 0L,
                            ifelse(is.na(negative_earliest_date), 1L,
                                   ifelse(!is.na(positive_earliest_date) & !is.na(negative_earliest_date) & negative_latest_date<positive_earliest_date, 1L, NA)))) %>%
  analysis$cached("final_fh_24", unique_indexes="patid")

cohort <- cohort %>%
  left_join((final_fh %>% select(patid, fh_diabetes)), by="patid") %>%
  analysis$cached("cohort_24_interim_12", unique_indexes="patid")


############################################################################################

# Add earliest insulin and OHA (pre-index date)

## Earliest insulin per patient

earliest_insulin <- clean_insulin_prodcodes %>%
  filter(date<=index_date) %>%
  group_by(patid) %>%
  summarise(earliest_ins=min(date, na.rm=TRUE)) %>%
  ungroup() %>%
  analysis$cached("earliest_insulin_24", unique_indexes="patid")


## Earliest OHA per patient

earliest_oha <- clean_oha_prodcodes %>%
  filter(date<=index_date) %>%
  group_by(patid) %>%
  summarise(earliest_oha=min(date, na.rm=TRUE)) %>%
  ungroup() %>%
  analysis$cached("earliest_oha_24", unique_indexes="patid")


## Combine and add time to insulin within 3 years

cohort <- cohort %>%
  left_join(earliest_insulin, by="patid") %>%
  left_join(earliest_oha, by="patid") %>%
  analysis$cached("cohort_24", unique_indexes="patid")

cohort %>% count()
#1,278,300