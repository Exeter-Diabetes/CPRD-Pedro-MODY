# CPRD-Pedro-MODY
R scripts for MODY prediction in CPRD


## MODY T1D model

Bayesian hierarchical model combining case-control data with UNITED (population representative data). It uses a mixture approach splitting patients based on a latent variable *T*. If the patient is $C^- \cup A^+$ then $T = 1$. If the patient is $C^+ \cap A^-$ then $T = 0$. Patients with $T = 1$ have an informative prior Beta distribution, modelling the probability of MODY, with $mean = 0.029%$. Patients with $T = 0$ are modelled using shrinkage recalibration logistic model which scales the odds ratios and adjusts the intercept based on the general population data. In case *T* is missing, *T* is modelled using a logistic regression with the variables BMI, age of diagnosis, age of recruitment and parent history of diabetes (continuous variables are modelled with a restricted cubic spline, 3 knots calculated using UNITED dataset).


## MODY T2D model

Bayesian shrinkage recalibration logistic model which scales the odds ratios and adjusts the intercept based on the general population data.
