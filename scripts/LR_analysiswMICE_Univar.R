## Univariate log regression of depression score on region of chest pain
## Imputation done in separate code (see LR_analysiswMICE)

responses <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")
pooled_results <- list()

for (response in responses) {  
  fit_reg <- function(dat) {
    # Define survey design
    svy <- svydesign(
      id = ~SDMVPSU,
      weights = ~MEC15YR,
      strata = ~SDMVSTRA,
      data = dat,
      nest = TRUE
    )
    
    # Subset to exclude rows where all response variables are 0 or NA
    svy <- subset(svy, !apply(
      svy$variables[, responses], 
      1, 
      function(x) all(x == 0 | is.na(x))
    ))
    
    # Fit logistic regression model for the current response
    formula <- as.formula(paste(response, "~ DEPR_TOT"))
    mod <- svyglm(formula, design = svy, family = quasibinomial(link = "logit"))
    
    return(mod)
  }
  
  # Obtain all imputed datasets
  dat_mice <- complete(imputations, "all")
  
  # Apply the function to each imputed dataset
  mod_imputation <- lapply(dat_mice, fit_reg)
  
  # Pool the results across imputations
  mod_pooled <- pool(mod_imputation)
  
  # Store the pooled results for the current response
  pooled_results[[response]] <- summary(mod_pooled)
}

# Display pooled results
for (response in responses) {
  cat("Pooled Results for", response, ":\n")
  print(pooled_results[[response]])
  cat("\n")
}


















# Parallel backend setup
num_cores <- parallel::detectCores() - 1  # Use one less than the total cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)


# Parallelized analysis for DEPR_TOT as the sole predictor
simple_results <- foreach(imp_id = 1:imputations$m, .packages = c("survey")) %dopar% {
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
  
  # Initialize results for this imputation
  coef_list <- list()
  
  for (response in responses) {
    # Construct the formula
    formula <- as.formula(paste(response, "~ DEPR_TOT"))
    
    # Base model
    base_model <- svyglm(formula, design = rep_design, family = quasibinomial(link = "logit"))
    
    # Extract coefficient for DEPR_TOT
    base_coef <- coef(base_model)["DEPR_TOT"]
    
    # Bootstrap replicate analysis
    rep_coefficients <- numeric(ncol(weights(rep_design, type = "analysis")))
    
    for (i in seq_len(ncol(weights(rep_design, type = "analysis")))) {
      w <- weights(rep_design, type = "analysis")[, i]
      model <- svyglm(formula, design = rep_design, weights = w, family = quasibinomial(link = "logit"))
      rep_coefficients[i] <- coef(model)["DEPR_TOT"]
    }
    
    bootstrap_variance <- var(rep_coefficients)
    percentile_ci <- quantile(rep_coefficients, probs = c(0.025, 0.975))
    
    results <- data.frame(
      Response = response,
      Coefficient = "DEPR_TOT",
      Base_Coef = base_coef,
      Bootstrap_Var = bootstrap_variance,
      CI_Lower = percentile_ci[1],
      CI_Upper = percentile_ci[2]
    )
    coef_list[[response]] <- results
  }
  
  coef_list
}

# Pool results manually using Rubin's rules
pooled_simple_results <- list()

for (response in responses) {
  # Collect coefficients across imputations
  coefficients <- sapply(simple_results, function(x) x[[response]]$Base_Coef)
  bootstrap_vars <- sapply(simple_results, function(x) x[[response]]$Bootstrap_Var)
  
  # Mean coefficients (point estimates)
  mean_coef <- mean(coefficients)
  
  # Within-imputation variance
  within_var <- mean(bootstrap_vars)
  
  # Between-imputation variance
  between_var <- var(coefficients)
  
  # Total variance
  total_var <- within_var + (1 + 1 / length(coefficients)) * between_var
  pooled_se <- sqrt(total_var)
  
  # Confidence intervals
  ci_lower <- mean_coef - 1.96 * pooled_se
  ci_upper <- mean_coef + 1.96 * pooled_se
  
  # Convert to OR
  or <- exp(mean_coef)
  ci_lower_or <- exp(ci_lower)
  ci_upper_or <- exp(ci_upper)
  
  pooled_simple_results[[response]] <- data.frame(
    Coefficient = "DEPR_TOT",
    Pooled_OR = or,
    Pooled_OR_CI_Lower = ci_lower_or,
    Pooled_OR_CI_Upper = ci_upper_or,
    Pooled_SE = pooled_se
  )
}

# Display results
for (response in responses) {
  print(pooled_simple_results[[response]])
}
