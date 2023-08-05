############################################################################
############################################################################
############################################################################
#
#       This file contains the prediction functions for the MODY models 
#
############################################################################
############################################################################
############################################################################


# load libraries
require(tidyverse)
require(nimble)

## Prediction function for MODY T1D

### predict method for 'post' objects
predict.T1D <- function(object, newdata, parms, ...) {
  
  ## check input objects
  stopifnot(class(object) == "T1D")
  stopifnot(is.data.frame(newdata))
  stopifnot(c("pardm", "agerec", "hba1c", "agedx", "sex", "bmi", "T") %in% colnames(newdata))
  x <- select(newdata, one_of(c("pardm", "agerec", "hba1c", "agedx", "sex")))
  xT <- as.matrix(select(newdata, bmi, agedx, pardm, agerec))
  x_spline <- newdata %>%
    select(bmi, agedx, agerec) %>%
    mutate(bmi_spline = as.numeric(rms::rcs(bmi, parms = parms[,"bmi"])[, 2]),
           agedx_spline = as.numeric(rms::rcs(agedx, parms = parms[,"agedx"])[, 2]),
           agerec_spline = as.numeric(rms::rcs(agerec, parms = parms[,"agerec"])[, 2])) %>%
    select(-bmi, -agedx, -agerec) %>%
    set_names(c("bmi_spline", "agedx_spline", "agerec_spline")) %>%
    as.matrix()
  stopifnot(all(map_chr(x, class) == "numeric"))
  stopifnot(is.numeric(newdata$T) | all(is.na(newdata$T)))
  
  ## convert posterior samples to matrix
  post <- as.matrix(object$post)
  
  # ## set up data## set up data
  x <- as.matrix(x)
  x <- cbind(int = rep(1, nrow(x)), x)
  betas <- post[, match(c("beta0", paste0("beta[", 1:(ncol(x) - 1), "]")), colnames(post))]
  betas <- x %*% t(betas)
  
  ## do predictions ignoring tests
  preds <- t(t(betas) * post[, "gamma1"])
  preds <- t(t(preds) + post[, "gamma0"])
  preds <- exp(preds) / (1 + exp(preds))
  
  ## all possible combinations of tests
  known_T <- which(!is.na(newdata$T))
  known_nT <- which(is.na(newdata$T))
  
  ## adjust predictions for known T
  if (sum(known_T) > 0) {
    
    # collect values of T and predictions
    T <- newdata$T[!is.na(newdata$T)]
    preds_T <- preds[!is.na(newdata$T),]
    
    # add predictions weighted by variable T
    pT <- matrix(T, ncol = 1) %*% t(post[, "pMp_Cn_or_Ap"]) + matrix(replicate(ncol(preds_T), 1 - T), ncol = ncol(preds_T)) * preds_T
    
    # replace new probability calculations
    preds[!is.na(newdata$T)] <- pT
    
  }
  
  ## adjust predictions for not known T
  if (sum(known_nT) >0) {
    
    # collect values predictions
    preds_nT <- preds[is.na(newdata$T),]
    
    # calculate variable T
    x_all <- cbind(int = rep(1, nrow(x)), xT, x_spline)
    betas <- post[, match(c("beta_t0", paste0("beta_t[", 1:ncol(xT), "]"), paste0("beta_spline[", 1:ncol(x_spline), "]")), colnames(post))]
    betas <- x_all %*% t(betas)
    
    pT_result <- exp(betas) / (1 + exp(betas))
    pT_result <- pT_result[is.na(newdata$T),]
    
    # add predictions weighted by variable T
    if (length(known_nT) == 1) {
      pnT <- pT_result * t(replicate(1, post[, "pMp_Cn_or_Ap"])) + (1 - pT_result)  * preds_nT
    } else {
      pnT <- pT_result * t(replicate(nrow(pT_result), post[, "pMp_Cn_or_Ap"])) + (1 - pT_result)  * preds_nT
    }
    
    # replace new probability calculations
    preds[is.na(newdata$T)] <- pnT
    
  }
  
  ## return posterior predictive samples
  preds <- t(preds)
  colnames(preds) <- NULL
  preds
}


## Prediction function for MODY T2D

### predict method for 'post' objects
predict.T2D <- function(object, newdata, ...) {
  
  ## check input objects
  stopifnot(class(object) == "T2D")
  stopifnot(is.data.frame(newdata))
  stopifnot(c("agedx", "bmi", "hba1c", "pardm", "agerec", "insoroha", "sex") %in% colnames(newdata))
  x <- select(newdata, one_of(c("agedx", "bmi", "hba1c", "pardm", "agerec", "insoroha", "sex")))
  stopifnot(all(map_chr(x, class) == "numeric"))
  
  ## convert posterior samples to matrix
  post <- as.matrix(object$post)
  
  ## set up data
  x <- as.matrix(x)
  x <- cbind(int = rep(1, nrow(x)), x)
  betas <- post[, match(c("beta0", paste0("beta[", 1:(ncol(x) - 1), "]")), colnames(post))]
  betas <- x %*% t(betas)
  
  ## do predictions ignoring tests
  preds <- t(t(betas) * post[, "gamma1"])
  preds <- t(t(preds) + post[, "gamma0"])
  preds <- exp(preds) / (1 + exp(preds))
  
  ## return posterior predictive samples
  preds <- t(preds)
  colnames(preds) <- NULL
  preds
}

