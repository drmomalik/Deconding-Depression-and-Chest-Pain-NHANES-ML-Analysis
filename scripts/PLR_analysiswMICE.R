# Load required libraries
library(mice)
library(glmnet)
library(pROC)
library(doParallel)
library(foreach)
library(tidyverse)
library(survey)


# Parallel backend setup
num_cores <- parallel::detectCores() - 1  # Use one less than the total cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Subset variables to improve computational efficiency
responses <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")
predictors <- c("DEPR_BIN", "DEPR_LVL", "AGE_BIN", "RIAGENDR", "RIDRETH1", "SMQ020",
                "SMQ040", "HIQ011", "BPQ020", "BMXBMI", "BMI_LVL", "BPQ080",
                "ALQ130", "MEDDEP", "DMDBORNT", "PAQMV", "CADTOT", "DIDTOT",
                "DUQTOT", "INC_BIN", "INC3", "DMDEDUC2", "CDQ008", "CDQ010")
wd_subset <- wd[, c(responses, predictors,"SDDSRVYR", "SDMVPSU", "SDMVSTRA", "MEC15YR", "SEQN")]

# Initialize the predictor matrix for MICE
predictorMatrix <- mice::make.predictorMatrix(wd_subset)
predictorMatrix[, !colnames(predictorMatrix) %in% predictors] <- 0
predictorMatrix[!rownames(predictorMatrix) %in% predictors, ] <- 0

# Run MICE
imputations <- mice(wd_subset, seed = 123, maxit = 1, m = 1, predictorMatrix = predictorMatrix, printFlag = TRUE)

# Parallelized analysis for each imputed dataset
all_results <- foreach(imp_id = 1:imputations$m, .packages = c("glmnet", "survey")) %dopar% {
  # Extract the imputed dataset
  imputation <- mice::complete(imputations, imp_id)
  
  # Survey design setup
  design <- svydesign(
    id = ~SDMVPSU,
    weights = ~MEC15YR,
    strata = ~SDMVSTRA,
    data = imputation,
    nest = TRUE
  )
  set.seed(123)
  rep_design <- as.svrepdesign(design, type = "bootstrap", replicates = 500)
  rep_design <- subset(rep_design, !apply(
    rep_design$variables[, responses], 1,
    function(x) all(x == 0 | is.na(x))
  ))
  
  # Stabilize weights
  repweights <- weights(rep_design, type = "analysis") / mean(weights(rep_design, type = "analysis"))
  pweights <- weights(rep_design, type = "sampling") / mean(weights(rep_design, type = "sampling"))
  
  # Transform variables
  rep_design$variables$BMXBMI <- (rep_design$variables$BMXBMI)^2
  
  # Create the model matrix
  X <- model.matrix(
    ~ DEPR_BIN + AGE_BIN + RIAGENDR + RIDRETH1 + SMQ040 + SMQ020 +
      HIQ011 + BPQ020 + BMI_LVL + BPQ080 + ALQ130 + MEDDEP + DMDBORNT + 
      PAQMV + CADTOT + DIDTOT + DUQTOT + INC_BIN +
      SMQ040 * ALQ130 * DUQTOT + PAQMV * BMI_LVL + DEPR_BIN * MEDDEP + DEPR_BIN * CADTOT,
    data = rep_design$variables
  )
  
  # Initialize results for this imputation
  coef_list <- list()
  
  for (response in responses) {
    y <- rep_design$variables[[response]]
    
    # Create penalty factors: lower penalty for priority_vars
    priority_vars <- c("DEPR_BIN2", "DEPR_BIN2:MEDDEP1", "DEPR_BIN2:CADTOT1")
    penalty_factors <- ifelse(colnames(X) %in% priority_vars, 0, 1)
    
    base_model <- cv.glmnet(
      X, y, family = "binomial", weights = pweights, penalty.factor = penalty_factors,
      alpha = 0, nfolds = 10, type.measure = "auc", standardize = TRUE
    )
    
    # Calculate marginal effect and coefficients
    probs <- table(rep_design$variables$MEDDEP, rep_design$variables$CADTOT) / nrow(rep_design$variables)
    marginal_effect <- coef(base_model, s = "lambda.min")@x[2] + 
      coef(base_model, s = "lambda.min")@x[48] * probs["1", "0"] +
      coef(base_model, s = "lambda.min")@x[49] * probs["0", "1"] + 
      (coef(base_model, s = "lambda.min")@x[48] + coef(base_model, s = "lambda.min")@x[49]) * probs["1", "1"]
    
    base_coef <- c(marginal_effect, coef(base_model, s = "lambda.min")@x[48], coef(base_model, s = "lambda.min")@x[49])
    names(base_coef) <- c("DEPR_BIN", "DEPR_BIN:MEDDEP", "DEPR_BIN:CADTOT")
    
    # Bootstrap replicate analysis
    rep_coef_matrix <- matrix(NA, nrow = 3, ncol = ncol(repweights))
    rownames(rep_coef_matrix) <- c("DEPR_BIN", "DEPR_BIN:MEDDEP", "DEPR_BIN:CADTOT")
    
    for (i in seq_len(ncol(repweights))) {
      w <- repweights[, i]
      model <- cv.glmnet(
        X, y, family = "binomial", weights = w, penalty.factor = penalty_factors,
        alpha = 0, nfolds = 10, type.measure = "auc", standardize = TRUE
      )
      marginal_effect <- coef(model, s = "lambda.min")@x[2] + 
        coef(model, s = "lambda.min")@x[48] * probs["1", "0"] +
        coef(model, s = "lambda.min")@x[49] * probs["0", "1"] + 
        (coef(model, s = "lambda.min")@x[48] + coef(model, s = "lambda.min")@x[49]) * probs["1", "1"]
      rep_coef_matrix[, i] <- c(marginal_effect, coef(model, s = "lambda.min")@x[48], coef(model, s = "lambda.min")@x[49])
    }
    
    rep_coef_df <- as.data.frame(rep_coef_matrix)
    colnames(rep_coef_df) <- paste0("Replicate_", seq_len(ncol(repweights)))
    bootstrap_variance <- apply(rep_coef_df, 1, var)
    percentile_ci <- t(apply(rep_coef_df, 1, quantile, probs = c(0.025, 0.975)))
    
    results <- data.frame(
      Response = response,
      Coefficient = rownames(rep_coef_df),
      Base_Coef = base_coef,
      Bootstrap_Var = bootstrap_variance,
      CI_Lower = percentile_ci[, 1],
      CI_Upper = percentile_ci[, 2]
    )
    coef_list[[response]] <- results
  }
  
  coef_list
}

# Stop the parallel cluster
stopCluster(cl)

# After running the foreach loop and collecting results in `all_results` (a list of lists)
# We need to pool the results manually.

# Initialize a list to store pooled results
pooled_results <- list()

# Loop over each response to apply Rubin's rules
for (response in responses) {
  # Extract coefficients and standard errors for this response across all imputations
  coefs <- sapply(all_results, function(x) x[[response]]$Base_Coef)
  se <- sapply(all_results, function(x) x[[response]]$Bootstrap_Var)
  
  # Calculate the mean across imputations (point estimate)
  mean_coef <- rowMeans(coefs, na.rm = TRUE)
  
  # Calculate within-imputation variance (average of the individual variances)
  within_var <- rowMeans(se, na.rm = TRUE)
  
  # Calculate between-imputation variance (variance of the point estimates across imputations)
  between_var <- apply(coefs, 1, var, na.rm = TRUE)
  
  # Calculate total variance (within + between)
  total_var <- within_var + (1 + 1/length(all_results)) * between_var
  
  # Calculate the pooled standard error
  pooled_se <- sqrt(total_var)
  
  # Calculate the pooled coefficients (mean of the coefficients)
  pooled_coef <- mean_coef
  
  # Compute 95% confidence intervals (using the pooled SE and normal approximation)
  ci_lower <- pooled_coef - 1.96 * pooled_se
  ci_upper <- pooled_coef + 1.96 * pooled_se
  
  # Store the results in a data frame
  pooled_results[[response]] <- data.frame(
    Coefficient = rownames(coefs),
    Pooled_Coef = pooled_coef,
    Pooled_SE = pooled_se,
    CI_Lower = ci_lower,
    CI_Upper = ci_upper
  )
}

# Display pooled results
for (response in responses) {
  print(pooled_results[[response]])
}
