############################################################################
############################################################################
############################################################################
#
#       In this file we load CPRD predictions and make a file that
#         can be passed around with just predictions.
#
############################################################################
############################################################################
############################################################################

# load libraries
require(tidyverse)

# load dataset - modyt1d_cohort_local
load("/slade/CPRD_data/Katie Pedro MODY/pedro_mody_cohort_2024_v2.Rda")

## remove those with missing bmi, hba1c
modyt1d_cohort_local_clean <- pedro_mody_cohort_local %>%
  drop_na(bmi, hba1c)

# load calculator predictions
T1D_predictions_no_T <- readRDS("CPRD_MODY/Patient Predictions/T1D_predictions_no_T.rds")
T2D_predictions <- readRDS("CPRD_MODY/Patient Predictions/T2D_predictions.rds")


#############################################
## Predictions from old calculator

### predict method for 'post' objects
predict.old_calculator_T1D <- function(object, newdata, ...) {
  
  ## check input objects
  stopifnot(class(object) == "old_calculator_T1D")
  stopifnot(is.data.frame(newdata))
  stopifnot(c("pardm", "agerec", "hba1c", "agedx", "sex") %in% colnames(newdata))
  x <- select(newdata, one_of(c("pardm", "agerec", "hba1c", "agedx", "sex")))
  stopifnot(all(map_chr(x, class) == "numeric"))
  
  ## convert posterior samples to matrix
  post <- as.matrix(do.call(rbind, object$post))
  
  ## set up data
  x <- as.matrix(x)
  x <- cbind(int = rep(1, nrow(x)), x)
  betas <- post[, match(c("beta0", paste0("beta[", 1:(ncol(x) - 1), "]")), colnames(post))]
  betas <- x %*% t(betas)
  
  ## do predictions ignoring tests
  preds <- exp(betas) / (1 + exp(betas))
  
  ## return posterior predictive samples
  preds <- t(preds)
  colnames(preds) <- NULL
  preds
}

### predict method for 'post' objects
predict.old_calculator_T2D <- function(object, newdata, ...) {
  
  ## check input objects
  stopifnot(class(object) == "old_calculator_T2D")
  stopifnot(is.data.frame(newdata))
  stopifnot(c("agedx", "bmi", "hba1c", "pardm", "agerec", "insoroha", "sex") %in% colnames(newdata))
  x <- select(newdata, one_of(c("agedx", "bmi", "hba1c", "pardm", "agerec", "insoroha", "sex")))
  stopifnot(all(map_chr(x, class) == "numeric"))
  
  ## convert posterior samples to matrix
  post <- as.matrix(do.call(rbind, object$post))
  
  ## set up data
  x <- as.matrix(x)
  x <- cbind(int = rep(1, nrow(x)), x)
  betas <- post[, match(c("beta0", paste0("beta[", 1:(ncol(x) - 1), "]")), colnames(post))]
  betas <- x %*% t(betas)
  
  ## do predictions ignoring tests
  preds <- exp(betas) / (1 + exp(betas))
  
  ## return posterior predictive samples
  preds <- t(preds)
  colnames(preds) <- NULL
  preds
}


## risk conversion (original method)
convert <- tibble(
  threshold = round(seq(0, 0.9, by = 0.1) * 100, 0), 
  PPVT1 = c(0.007, 0.019, 0.026, 0.040, 0.049, 0.064, 0.072, 0.082, 0.126, 0.494),
  PPVT2 = c(0.046, 0.151, 0.210, 0.244, 0.329, 0.358, 0.455, 0.580, 0.624, 0.755)
)

## Old calculator
posteriors_samples_old_T1D <- readRDS("CPRD_MODY/model_posteriors/type_1_old_model_posteriors.rds")
### create object to use for prediction
posteriors_samples_old_T1D <- list(post = posteriors_samples_old_T1D$samples)
class(posteriors_samples_old_T1D) <- "old_calculator_T1D"

posteriors_samples_old_T2D <- readRDS("CPRD_MODY/model_posteriors/type_2_old_model_posteriors.rds")
### create object to use for prediction
posteriors_samples_old_T2D <- list(post = posteriors_samples_old_T2D$samples)
class(posteriors_samples_old_T2D) <- "old_calculator_T2D"


### Type 1 old model

## make predictions for T1D MODY (missingness is in pardm)
final_T1D_predictions <- data.frame(
  patid = modyt1d_cohort_local_clean$patid
)

interim <- as_tibble(as.matrix(select(modyt1d_cohort_local_clean, pardm, agerec, hba1c, agedx, sex))) %>%
  mutate(pardm = ifelse(is.na(pardm), 0, pardm))

T1D_predictions_old <- predict(posteriors_samples_old_T1D, interim) %>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x))
  }) %>%
  bind_rows()

for (i in 1:nrow(T1D_predictions_old)) {
  prob <- T1D_predictions_old[i,]
  prob <- 10 * ((round(prob  * 100, 0)) %/% 10)
  if (prob == 100) {
    T1D_predictions_old[i,] <- convert$PPVT1[10]
  } else {
    T1D_predictions_old[i,] <- convert$PPVT1[convert$threshold == prob]
  }
}

final_T1D_predictions <- final_T1D_predictions %>%
  left_join(
    T1D_predictions_old %>%
      cbind(
        patid = modyt1d_cohort_local_clean$patid
      ) %>%
      set_names(c("mean_pardm_0", "patid"))
  )



interim <- as_tibble(as.matrix(select(modyt1d_cohort_local_clean, pardm, agerec, hba1c, agedx, sex))) %>%
  mutate(pardm = ifelse(is.na(pardm), 1, pardm))

T1D_predictions_old <- predict(posteriors_samples_old_T1D, interim) %>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x))
  }) %>%
  bind_rows()

for (i in 1:nrow(T1D_predictions_old)) {
  prob <- T1D_predictions_old[i,]
  prob <- 10 * ((round(prob  * 100, 0)) %/% 10)
  if (prob == 100) {
    T1D_predictions_old[i,] <- convert$PPVT1[10]
  } else {
    T1D_predictions_old[i,] <- convert$PPVT1[convert$threshold == prob]
  }
}

final_T1D_predictions <- final_T1D_predictions %>%
  left_join(
    T1D_predictions_old %>%
      cbind(
        patid = modyt1d_cohort_local_clean$patid
      ) %>%
      set_names(c("mean_pardm_1", "patid"))
  )


T1D_predictions_old <- final_T1D_predictions


### Type 2 old model

## make predictions for T1D MODY (missingness is in pardm)
final_T2D_predictions <- data.frame(
  patid = modyt1d_cohort_local_clean$patid
)

interim <- as_tibble(as.matrix(select(modyt1d_cohort_local_clean, agedx, bmi, hba1c, pardm, agerec, insoroha, sex))) %>%
  mutate(pardm = ifelse(is.na(pardm), 0, pardm))

T2D_predictions_old <- predict(posteriors_samples_old_T2D, interim) %>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x))
  }) %>%
  bind_rows()

for (i in 1:nrow(T2D_predictions_old)) {
  prob <- T2D_predictions_old[i,]
  prob <- 10 * ((round(prob  * 100, 0)) %/% 10)
  if (prob == 100) {
    T2D_predictions_old[i,] <- convert$PPVT2[10]
  } else {
    T2D_predictions_old[i,] <- convert$PPVT2[convert$threshold == prob]
  }
}

final_T2D_predictions <- final_T2D_predictions %>%
  left_join(
    T2D_predictions_old %>%
      cbind(
        patid = modyt1d_cohort_local_clean$patid
      ) %>%
      set_names(c("mean_pardm_0", "patid"))
  )



interim <- as_tibble(as.matrix(select(modyt1d_cohort_local_clean, agedx, bmi, hba1c, pardm, agerec, insoroha, sex))) %>%
  mutate(pardm = ifelse(is.na(pardm), 1, pardm))

T2D_predictions_old <- predict(posteriors_samples_old_T2D, interim) %>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x))
  }) %>%
  bind_rows()

for (i in 1:nrow(T2D_predictions_old)) {
  prob <- T2D_predictions_old[i,]
  prob <- 10 * ((round(prob  * 100, 0)) %/% 10)
  if (prob == 100) {
    T2D_predictions_old[i,] <- convert$PPVT2[10]
  } else {
    T2D_predictions_old[i,] <- convert$PPVT2[convert$threshold == prob]
  }
}

final_T2D_predictions <- final_T2D_predictions %>%
  left_join(
    T2D_predictions_old %>%
      cbind(
        patid = modyt1d_cohort_local_clean$patid
      ) %>%
      set_names(c("mean_pardm_1", "patid"))
  )


T2D_predictions_old <- final_T2D_predictions



#############################################

## Creating tables needed

## Take the original dataset and keep patid and calculator needed + agedx

mody_prob_CPRD <- modyt1d_cohort_local_clean %>%
  select(patid, diabetes_type, agedx) %>%
  # join old calculator probability T1D
  left_join(
    T1D_predictions_old %>%
      # create final probability = average of pardm = 1 and pardm = 0
      mutate(mean_T1D_final = mean_pardm_0) %>%
      select(patid, mean_T1D_final),
    by = c("patid")
  ) %>%
  # join old calculator probability T2D
  left_join(
    T2D_predictions_old %>%
      # create final probability = average of pardm = 1 and pardm = 0
      mutate(mean_T2D_final = mean_pardm_0) %>%
      select(patid, mean_T2D_final),
    by = c("patid")
  ) %>%
  # select the right calculator for each patient
  mutate(old_prob = ifelse(grepl('type 1',diabetes_type), mean_T1D_final, mean_T2D_final)) %>%
  # keep only important columns
  select(patid, diabetes_type, agedx, old_prob) %>%
  # join new calculator probability T1D
  left_join(
    T1D_predictions_no_T %>%
      # create final probability = average of pardm = 1 and pardm = 0
      mutate(mean_T1D_final = mean_pardm_0) %>%
      select(patid, mean_T1D_final),
    by = c("patid")
  ) %>%
  # join new calculator probability T2D
  left_join(
    T2D_predictions %>%
      # create final probability = average of pardm = 1 and pardm = 0
      mutate(mean_T2D_final = mean_pardm_0) %>%
      select(patid, mean_T2D_final),
    by = c("patid")
  ) %>%
  # select the right calculator for each patient
  mutate(new_prob = ifelse(grepl('type 1',diabetes_type), mean_T1D_final, mean_T2D_final)) %>%
  # keep only important columns
  select(patid, diabetes_type, agedx, old_prob, new_prob) %>%
  # simpler version of which calculator
  mutate(which_eq = ifelse(grepl('type 1',diabetes_type), "t1d", "t2d")) %>%
  select(-diabetes_type)
  

## under 35s table
under_35 <- mody_prob_CPRD %>%
  select(which_eq, old_prob, new_prob)

write.table(under_35, file = "CPRD_MODY/Patient Predictions/mody_probabilities_under_35s.txt", sep = "\t",
            row.names = FALSE)

## under 30s table
under_30 <- mody_prob_CPRD %>%
  filter(agedx < 30) %>%
  select(which_eq, old_prob, new_prob)

write.table(under_30, file = "CPRD_MODY/Patient Predictions/mody_probabilities_under_30s.txt", sep = "\t",
            row.names = FALSE)















