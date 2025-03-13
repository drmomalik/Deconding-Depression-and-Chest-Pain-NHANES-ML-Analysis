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
imputations <- mice(wd_subset, seed = 123, maxit = 5, m = 5, predictorMatrix = predictorMatrix, printFlag = TRUE)

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
  
  # Initialize results for this imputation
  coef_list <- list()
  
  for (response in responses) {
    pred <- paste(
      "AGE_BIN + RIAGENDR + RIDRETH1 +",
      "HIQ011 + I(BMXBMI^2) + INC3 + DMDEDUC2 +",
      "PAQMV + I(ALQ130^2) + DIDTOT + DUQTOT + SMQ040 * SMQ020 +",
      "DMDBORNT + BPQ020 + BPQ080 + SDDSRVYR +",
      "DEPR_BIN * MEDDEP + DEPR_BIN * CADTOT +",
      "CDQ008 + CDQ010"
    )
    
    # Construct the formula
    formula <- as.formula(paste(response, "~", pred))
    
    # Base model
    base_model <- glm(formula, 
                      data = rep_design$variables, weights = pweights, 
                      family = quasibinomial(link = "logit"))
    
    # Extract probabilities for MEDDEP and CADTOT
    probs <- table(rep_design$variables$MEDDEP, rep_design$variables$CADTOT) / nrow(rep_design$variables)
    
    # Marginal effects for DEPR_BIN
    base_marginals <- with(as.list(coef(base_model)), {
      main_effect <- DEPR_BIN2
      meddep_interaction <- `DEPR_BIN2:MEDDEP1`
      cadtot_interaction <- `DEPR_BIN2:CADTOT1`
      main_effect + 
        meddep_interaction * probs["1", "0"] +
        cadtot_interaction * probs["0", "1"] +
        (meddep_interaction + cadtot_interaction) * probs["1", "1"]
    })
    
    # Bootstrap replicate analysis
    rep_marginals <- numeric(ncol(repweights))
    
    for (i in seq_len(ncol(repweights))) {
      w <- repweights[, i]
      model <- glm(formula,
                   data = rep_design$variables, weights = w, 
                   family = quasibinomial(link = "logit"))
      rep_marginals[i] <- with(as.list(coef(model)), {
        main_effect <- DEPR_BIN2
        meddep_interaction <- `DEPR_BIN2:MEDDEP1`
        cadtot_interaction <- `DEPR_BIN2:CADTOT1`
        main_effect + 
          meddep_interaction * probs["1", "0"] +
          cadtot_interaction * probs["0", "1"] +
          (meddep_interaction + cadtot_interaction) * probs["1", "1"]
      })
    }
    
    bootstrap_variance <- var(rep_marginals)
    percentile_ci <- quantile(rep_marginals, probs = c(0.025, 0.975))
    
    results <- data.frame(
      Response = response,
      Coefficient = "DEPR_BIN",
      Base_Marginal = base_marginals,
      Bootstrap_Var = bootstrap_variance,
      CI_Lower = percentile_ci[1],
      CI_Upper = percentile_ci[2]
    )
    coef_list[[response]] <- results
  }
  
  coef_list
}

# Stop the parallel cluster
stopCluster(cl)

# Pool results manually using Rubin's rules
pooled_results <- list()

for (response in responses) {
  # Collect marginal effects across imputations
  marginal_effects <- sapply(all_results, function(x) x[[response]]$Base_Marginal)
  bootstrap_vars <- sapply(all_results, function(x) x[[response]]$Bootstrap_Var)
  
  # Mean coefficients (point estimates)
  mean_coef <- mean(marginal_effects)
  
  # Within-imputation variance
  within_var <- mean(bootstrap_vars)
  
  # Between-imputation variance
  between_var <- var(marginal_effects)
  
  # Total variance
  total_var <- within_var + (1 + 1 / length(marginal_effects)) * between_var
  pooled_se <- sqrt(total_var)
  
  # Confidence intervals
  ci_lower <- mean_coef - 1.96 * pooled_se
  ci_upper <- mean_coef + 1.96 * pooled_se
  
  pooled_results[[response]] <- data.frame(
    Coefficient = "DEPR_BIN",
    Pooled_Marginal = mean_coef,
    Pooled_SE = pooled_se,
    CI_Lower = ci_lower,
    CI_Upper = ci_upper
  )
}

# Display results
for (response in responses) {
  print(pooled_results[[response]])
}

