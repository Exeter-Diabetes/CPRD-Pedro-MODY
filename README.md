# CPRD-Pedro-MODY

R scripts for MODY prediction in CPRD

## MODY T1D model

Bayesian hierarchical model combining case-control data with UNITED (population representative data). It uses a mixture approach splitting patients based on a latent variable *T*. If the patient is $C^- \cup A^+$ then $T = 1$. If the patient is $C^+ \cap A^-$ then $T = 0$. Patients with $T = 1$ have an informative prior Beta distribution, modelling the probability of MODY, with $mean = 0.029%$. Patients with $T = 0$ are modelled using shrinkage recalibration logistic model which scales the odds ratios and adjusts the intercept based on the general population data. In case *T* is missing, *T* is modelled using a logistic regression with the variables BMI, age of diagnosis, age of recruitment and parent history of diabetes (continuous variables are modelled with a restricted cubic spline, 3 knots calculated using UNITED dataset).

### How to make predictions

The function used for predictions is available at `00.prediction_functions.R`. In order to use it, you must load the function into the environment.

```         
source("00.prediction_functions.R")
```

You also need to load the Bayesian model posteriors (model parameters) and format them the right way.

```         
rcs_parms <- readRDS("model_posteriors/rcs_parms.rds")
posterior_samples_T1D <- readRDS("model_posteriors/type_1_model_posteriors.rds")
posterior_samples_T1D_obj <- list(post = posterior_samples_T1D$samples)
class(posterior_samples_T1D_obj) <- "T1D"
```

In order to make predictions, you need the following variables: 
- `pardm`: Parent history of Diabetes, with history represented by `1` and no history represented by `0`
- `agerec`: Age at recruitment
- `hba1c`: HbA\_{1c}
- `agedx`: Age at diagnosis 
- `sex`: Sex, with Male represented by a `0` and Female represented by a `1`
- `bmi`: BMI
- `T`: C-peptide and autoantibody testing, with C+ and A- represented by `0`, C- or A+ represented by `1` and any other combination represented by `NA`.

The dataset has to be formatted as a tibble.

```
predictions_x <- as_tibble(as.matrix(select(patients_dataset, pardm, agerec, hba1c, agedx, sex, bmi, T)))
```

Last thing to do is the predictions themselves. For that, use

```
predictions_T1D <- predict(posterior_samples_T1D_obj, predictions_x, rcs_parms) %>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
  }) %>%
  bind_rows()
```

this will make the predictions and calculate the mean (lower and upper credible intervals at 2.5% and 97.5%) probability of having a MODY gene. If you are not looking at the uncertainty, only use the mean prediction.

## MODY T2D model

Bayesian shrinkage recalibration logistic model which scales the odds ratios and adjusts the intercept based on the general population data.

### How to make predictions

The function used for predictions is available at `00.prediction_functions.R`. In order to use it, you must load the function into the environment.

```         
source("00.prediction_functions.R")
```

You also need to load the Bayesian model posteriors (model parameters).

```         
posterior_samples_T2D <- readRDS("model_posteriors/type_2_model_posteriors.rds")
posterior_samples_T2D_obj <- list(post = posterior_samples_T2D$samples)
class(posterior_samples_T2D_obj) <- "T2D"
```

In order to make predictions, you need the following variables: 
- `pardm`: Parent history of Diabetes, with history represented by `1` and no history represented by `0`
- `agerec`: Age at recruitment
- `hba1c`: HbA\_{1c}
- `agedx`: Age at diagnosis 
- `sex`: Sex, with Male represented by a `0` and Female represented by a `1`
- `bmi`: BMI
- `insoroha`: Patient currently on insulin or tables, with TRUE represented by `1` and FALSE represented by `0`

The dataset has to be formatted as a tibble.

```
predictions_x <- as_tibble(as.matrix(select(patients_dataset, pardm, agerec, hba1c, agedx, sex, bmi, insoroha)))
```

Last thing to do is the predictions themselves. For that, use

```
predictions_T2D <- predict(posterior_samples_T2D_obj, predictions_x) %>%
  apply(., 2, function(x) {
    data.frame(prob = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
  }) %>%
  bind_rows()
```

this will make the predictions and calculate the mean (lower and upper credible intervals at 2.5% and 97.5%) probability of having a MODY gene. If you are not looking at the uncertainty, only use the mean prediction.
