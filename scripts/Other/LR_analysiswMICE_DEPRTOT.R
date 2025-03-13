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
predictors <- c("DEPR_TOT", "DEPR_BIN", "DEPR_LVL", "AGE_BIN", "RIAGENDR", "RIDRETH1", "SMQ020",
                "SMQ040", "HIQ011", "BPQ020", "BMXBMI", "BMI_LVL", "BPQ080",
                "ALQ130", "MEDDEP", "DMDBORNT", "PAQMV", "CADTOT", "DIDTOT",
                "DUQTOT", "INC_BIN", "INC3", "DMDEDUC2", "CDQ008", "CDQ010", "TOT_REG")
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
  rep_design <- as.svrepdesign(design, type = "bootstrap", replicates = 1000)
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
      "DEPR_TOT * MEDDEP + DEPR_TOT * CADTOT"
    )
    
    # Construct the formula
    formula <- as.formula(paste(response, "~", pred))
    
    # Base model
    base_model <- glm(formula, 
                      data = rep_design$variables, weights = pweights, 
                      family = quasibinomial(link = "logit"))
    
    # Extract probabilities for MEDDEP and CADTOT
    probs <- table(rep_design$variables$MEDDEP, rep_design$variables$CADTOT) / nrow(rep_design$variables)
    
    # Marginal effects for DEPR_TOT
    marginal_effect <- with(as.list(coef(base_model)), {
      main_effect <- DEPR_TOT
      meddep_interaction <- `DEPR_TOT:MEDDEP1`
      cadtot_interaction <- `DEPR_TOT:CADTOT1`
      main_effect + 
        meddep_interaction * probs["1", "0"] +
        cadtot_interaction * probs["0", "1"] +
        (meddep_interaction + cadtot_interaction) * probs["1", "1"]
    })
    
    # Bootstrap replicate analysis
    rep_marginals <- rep_main <- rep_meddep <- rep_cadtot <- numeric(ncol(repweights))
    
    for (i in seq_len(ncol(repweights))) {
      w <- repweights[, i]
      model <- glm(formula,
                   data = rep_design$variables, weights = w, 
                   family = quasibinomial(link = "logit"))
      rep_marginals[i] <- with(as.list(coef(model)), {
        main_effect <- DEPR_TOT
        meddep_interaction <- `DEPR_TOT:MEDDEP1`
        cadtot_interaction <- `DEPR_TOT:CADTOT1`
        main_effect + 
          meddep_interaction * probs["1", "0"] +
          cadtot_interaction * probs["0", "1"] +
          (meddep_interaction + cadtot_interaction) * probs["1", "1"]
      })
      rep_main[i] <- coef(model)["DEPR_TOT"]
      rep_meddep[i] <- coef(model)["DEPR_TOT:MEDDEP1"]
      rep_cadtot[i] <- coef(model)["DEPR_TOT:CADTOT1"]
    }
    
    # Calculate variance and CIs for marginal effect and coefficients
    bootstrap_variance <- c(
      Marginal = var(rep_marginals),
      DEPR_TOT = var(rep_main),
      Interaction_MEDDEP = var(rep_meddep),
      Interaction_CADTOT = var(rep_cadtot)
    )
    percentile_ci <- rbind(
      Marginal = quantile(rep_marginals, probs = c(0.025, 0.975)),
      DEPR_TOT = quantile(rep_main, probs = c(0.025, 0.975)),
      Interaction_MEDDEP = quantile(rep_meddep, probs = c(0.025, 0.975)),
      Interaction_CADTOT = quantile(rep_cadtot, probs = c(0.025, 0.975))
    )
    
    results <- data.frame(
      Response = response,
      Coefficient = c("Marginal_DEPR_TOT", "DEPR_TOT", "DEPR_TOT:MEDDEP", "DEPR_TOT:CADTOT"),
      Estimate = c(marginal_effect, coef(base_model)["DEPR_TOT"], coef(base_model)["DEPR_TOT:MEDDEP1"], coef(base_model)["DEPR_TOT:CADTOT1"]),
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

# Pool results manually using Rubin's rules
pooled_results <- list()

for (response in responses) {
  # Collect marginal effects and coefficients across imputations
  estimates <- sapply(all_results, function(x) x[[response]]$Estimate)
  bootstrap_vars <- sapply(all_results, function(x) x[[response]]$Bootstrap_Var)
  
  pooled_data <- t(apply(estimates, 1, function(coef_set) {
    mean_coef <- mean(coef_set)
    within_var <- mean(bootstrap_vars)
    between_var <- var(coef_set)
    total_var <- within_var + (1 + 1 / length(coef_set)) * between_var
    pooled_se <- sqrt(total_var)
    ci_lower <- mean_coef - 1.96 * pooled_se
    ci_upper <- mean_coef + 1.96 * pooled_se
    c(Mean = mean_coef, SE = pooled_se, CI_Lower = ci_lower, CI_Upper = ci_upper)
  }))
  
  pooled_results[[response]] <- as.data.frame(pooled_data)
  rownames(pooled_results[[response]]) <- c("Marginal_DEPR_TOT", "DEPR_TOT", "DEPR_TOT:MEDDEP", "DEPR_TOT:CADTOT")
}

# Display results
for (response in responses) {
  print(pooled_results[[response]])
}


# Exponentiate and flag significant rows
exp_pooled_results <- lapply(pooled_results, function(df) {
  df <- df %>%
    mutate(
      Mean = exp(Mean),
      CI_Lower = exp(CI_Lower),
      CI_Upper = exp(CI_Upper),
      Significant = (CI_Lower > 1 | CI_Upper < 1) # Flag rows where CI excludes 1
    )
  return(df)
})

# Function to print rows with highlighting for significant values
print_highlighted <- function(df) {
  apply(df, 1, function(row) {
    if (as.logical(row["Significant"])) {
      cat("\033[1m", paste(row, collapse = " | "), "\033[0m\n") # Bold the row
    } else {
      cat(paste(row, collapse = " | "), "\n")
    }
  })
}

# Display the results with highlighting
for (response in responses) {
  cat("Response:", response, "\n")
  print_highlighted(exp_pooled_results[[response]])
  cat("\n")
}