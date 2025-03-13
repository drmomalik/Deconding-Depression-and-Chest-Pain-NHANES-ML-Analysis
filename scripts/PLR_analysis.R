

### Penalized Log Reg
# Load glmnet library
library(glmnet)
library(pROC)

responses <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")

imputation <- complete(imputations, 1)

design <- svydesign(
  id = ~SDMVPSU, 
  weights = ~MEC15YR, 
  strata = ~SDMVSTRA, 
  data = imputation,
  nest = TRUE)

set.seed(123)  
rep_design <- as.svrepdesign(design, type = "bootstrap", replicates=100)
rep_design <- subset(rep_design, !apply(
  rep_design$variables[, c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")],
  1, 
  function(x) all(x == 0 | is.na(x))  # Keep rows where not all values of CDQ009A-CDQ009H are 0 or NA
))

# Create weights and stabilize 
repweights <- weights(rep_design, type = "analysis")
repweights <- repweights/mean(repweights)
pweights <- weights(rep_design, type = "sampling")
pweights <- pweights/mean(pweights)

## Variable transformation(s)
rep_design$variables$BMXBMI <- (rep_design$variables$BMXBMI)^2

# Prepare the data for glmnet
# Convert predictor variables to matrix form
X <- model.matrix(~ DEPR_BIN + AGE_BIN + RIAGENDR + RIDRETH1 + SMQ040 + SMQ020 +
                    HIQ011 + BPQ020 + BMI_LVL + BPQ080 + ALQ130 + MEDDEP + DMDBORNT + 
                    PAQMV + CADTOT + DIDTOT + DUQTOT + INC_BIN
                  + SMQ040*ALQ130*DUQTOT + PAQMV*BMI_LVL + DEPR_BIN*MEDDEP + DEPR_BIN*CADTOT
                  , data = rep_design$variables)

# Initialize a list for all responses

coef_list <- list()

for (response in responses) {

    # Convert the response variable to a vector
    y <- rep_design$variables[[response]]
    
    # Create penalty factors: lower penalty for DEPR_BIN2
    priority_vars <- c("DEPR_BIN2", "DEPR_BIN2:MEDDEP1", "DEPR_BIN2:CADTOT1")
    penalty_factors <- ifelse(colnames(X) %in% priority_vars, 0, 1)
    
    # Define the base model using glmnet
    base_model <- cv.glmnet(X, y, family = "binomial", weights = pweights, penalty.factor = penalty_factors,
                            alpha = 0, nfolds = 10, type.measure = "auc", standardize = TRUE)
    
    ## Find marginal effect of DEPR_BIN 
    probs <- table(rep_design$variables$MEDDEP, rep_design$variables$CADTOT) / nrow(rep_design$variables)
    marginal_effect <- coef(base_model, s="lambda.min")@x[2] + coef(base_model, s="lambda.min")@x[48]*probs["1","0"]
    + coef(base_model, s="lambda.min")@x[49]*probs["0","1"] + 
      (coef(base_model, s="lambda.min")@x[48]+coef(base_model, s="lambda.min")@x[49])*probs["1","1"]
    print(marginal_effect)
    
    
    # Extract coefficients for the base model and convert them to a vector
    base_coef <- c(marginal_effect, coef(base_model, s="lambda.min")@x[48], coef(base_model, s="lambda.min")@x[49])  # Use optimal lambda
    names(base_coef) <- c("DEPR_BIN", "DEPR_BIN:MEDDEP", "DEPR_BIN:CADTOT")
    
    # Initialize a matrix to store the coefficients
    rep_coef_matrix <- matrix(NA, nrow = 3, ncol = ncol(repweights))
    rownames(rep_coef_matrix) <- c("DEPR_BIN", "DEPR_BIN:MEDDEP", "DEPR_BIN:CADTOT")
    
    # Iterate over the replicate weights
    for (i in seq_len(ncol(repweights))) {
      # Extract the current replicate weight column
      w <- repweights[, i]
      
      # Fit the penalized logistic regression model for the current replicate weight
      model <- cv.glmnet(X, y, family = "binomial", weights = w, penalty.factor = penalty_factors,
                         alpha = 0, nfolds = 10, type.measure = "auc", standardize = TRUE)
      
      # Calculate the marginal effect for DEPR_BIN
      marginal_effect <- coef(model, s = "lambda.min")@x[2] + 
        coef(model, s = "lambda.min")@x[48] * probs["1", "0"] +
        coef(model, s = "lambda.min")@x[49] * probs["0", "1"] + 
        (coef(model, s = "lambda.min")@x[48] + coef(model, s = "lambda.min")@x[49]) * probs["1", "1"]
      
      # Extract the coefficients for the interactions
      depr_bin_meddep <- coef(model, s = "lambda.min")@x[48]
      depr_bin_cadtot <- coef(model, s = "lambda.min")@x[49]
      
      # Store the coefficients in the matrix
      rep_coef_matrix[, i] <- c(marginal_effect, depr_bin_meddep, depr_bin_cadtot)
    }
    
    # Convert the matrix to a dataframe for better readability
    rep_coef_df <- as.data.frame(rep_coef_matrix)
    colnames(rep_coef_df) <- paste0("Replicate_", seq_len(ncol(repweights)))
    
    # Display the resulting dataframe
    print(rep_coef_df)
    
    
    # Define a custom function for variance
    custom_var <- function(x) {
      mean_x <- mean(x)
      sum((x - mean_x)^2) / (length(x) - 1)
    }
    
    # Calculate bootstrap variance for each row
    bootstrap_variance <- apply(rep_coef_df, 1, custom_var)
    
    # Calculate the 95% percentile confidence intervals
    percentile_ci <- t(apply(rep_coef_df, 1, quantile, probs = c(0.025, 0.975)))
    
    # Combine the results into a dataframe
    results <- data.frame(
      Response = response,
      Coefficient = rownames(rep_coef_df),
      Base_Coef = base_coef,
      Bootstrap_Var = bootstrap_variance,
      CI_Lower = percentile_ci[, 1],
      CI_Upper = percentile_ci[, 2]
    )
    
    
    # Display the results
    print(results)
    
    coef_list[[response]] <- results    
}

#### 






