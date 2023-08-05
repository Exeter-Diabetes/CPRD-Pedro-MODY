############################################################################
############################################################################
############################################################################
#
#       In this file we load and format the CPRD dataset.
#
############################################################################
############################################################################
############################################################################

# load libraries
require(tidyverse)


#:---------------------------------------------------------

# load dataset - modyt1d_cohort_local
load("/slade/CPRD_data/mastermind_2022/Katie HF/pedro_modyt1d_cohort.Rda")
  
## remove those with missing bmi, hba1c
modyt1d_cohort_local_clean <- modyt1d_cohort_local %>%
  drop_na(bmi, hba1c)


#:---------------------------------------------------------

# load functions to make predictions from models
source("00.prediction_functions.R")

# load posteriors
rcs_parms <- readRDS("model_posteriors/rcs_parms.rds")
posterior_samples_T1D <- readRDS("model_posteriors/type_1_model_posteriors.rds")

### create object to use for prediction
posterior_samples_T1D_obj <- list(post = posterior_samples_T1D$samples)
class(posterior_samples_T1D_obj) <- "T1D"


#:---------------------------------------------------------

## make predictions for T1D MODY (missingness is in pardm)
final_T1D_predictions <- data.frame(
  patid = modyt1d_cohort_local_clean$patid
)



#:------------------
# If pardm is missing, set to 0
newdata_predictions <- modyt1d_cohort_local_clean %>% 
  mutate(pardm = ifelse(is.na(pardm), 0, pardm))

newdata_predictions_x <- as_tibble(as.matrix(select(newdata_predictions, pardm, agerec, hba1c, agedx, sex, bmi)))
newdata_predictions_x$T <- NA

predictions_T1D_pardm_0 <- predict(posterior_samples_T1D_obj, newdata_predictions_x, rcs_parms) %>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
  }) %>%
  bind_rows() %>%
  cbind(
    patid = newdata_predictions$patid
  )

final_T1D_predictions <- final_T1D_predictions %>%
  left_join(
    predictions_T1D_pardm_0 %>%
      set_names(c("mean_pardm_0", "lci_pardm_0", "uci_pardm_0", "patid"))
  )


#:------------------
# If pardm is missing, set to 1
newdata_predictions <- modyt1d_cohort_local_clean %>% 
  mutate(pardm = ifelse(is.na(pardm), 1, pardm))

newdata_predictions_x <- as_tibble(as.matrix(select(newdata_predictions, pardm, agerec, hba1c, agedx, sex, bmi)))
newdata_predictions_x$T <- NA

predictions_T1D_pardm_1 <- predict(posterior_samples_T1D_obj, newdata_predictions_x, rcs_parms) %>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
  }) %>%
  bind_rows() %>%
  cbind(
    patid = newdata_predictions$patid
  )

final_T1D_predictions <- final_T1D_predictions %>%
  left_join(
    predictions_T1D_pardm_1 %>%
      set_names(c("mean_pardm_1", "lci_pardm_1", "uci_pardm_1", "patid"))
  )


#:------------------
dir.create("Patient Predictions")
saveRDS(final_T1D_predictions, "Patient Predictions/T1D_predictions.rds")







