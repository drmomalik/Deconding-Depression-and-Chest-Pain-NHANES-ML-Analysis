# Survey design setup
design <- svydesign(
  id = ~SDMVPSU,
  weights = ~MEC15YR,
  strata = ~SDMVSTRA,
  data = wd_subset,
  nest = TRUE
)

survey_design <- subset(design, !apply(
  design$variables[, responses], 1,
  function(x) all(x == 0 | is.na(x))
))

survey_design$variables$SDMVSTRA <- NULL
survey_design$variables$SDMVPSU <- NULL
survey_design$variables$SEQN <- NULL

# missing Forest Imputation 
# 1. Subset the design into two groups based on CAD_TOT
survey_design_CAD1 <- subset(survey_design, CADTOT == 1)
survey_design_CAD0 <- subset(survey_design, CADTOT == 0)

survey_design_CAD1$variables$CADTOT <- NULL
survey_design_CAD0$variables$CADTOT <- NULL

# 2. Perform missing data imputation for each subset using missForest
# Convert the variables for imputation (subsetting just the variables)
subset_CAD1 <- survey_design_CAD1$variables
subset_CAD0 <- survey_design_CAD0$variables

# Impute missing values for each subset
imputed_CAD1 <- missForest(xmis = subset_CAD1, verbose = TRUE)
imputed_CAD0 <- missForest(xmis = subset_CAD0, verbose = TRUE)

# 3. Extract the imputed datasets
wd_subset_imputed_CAD1 <- imputed_CAD1$ximp  # Imputed data for CAD_TOT == 1
wd_subset_imputed_CAD0 <- imputed_CAD0$ximp  # Imputed data for CAD_TOT == 0

# 4. Add the imputed values back to the survey design objects
survey_design_CAD1$variables <- wd_subset_imputed_CAD1
survey_design_CAD0$variables <- wd_subset_imputed_CAD0

# Optionally, you can check the OOB errors for each imputation
print(imputed_CAD1$OOBerror)  # OOB error for CAD_TOT == 1
print(imputed_CAD0$OOBerror)  # OOB error for CAD_TOT == 0


train_survey_pcc_str <- function(survey_design, Y, label_order) {
  models <- list()
  pred <- c("DEPR_TOT","AGE_BIN", "RIAGENDR", "RIDRETH1", "SMQ020",
            "SMQ040", "HIQ011", "BPQ020", "BMXBMI", "BPQ080",
            "ALQ130", "MEDDEP", "DMDBORNT", "PAQMV", "DIDTOT",
            "DUQTOT", "INC3", "DMDEDUC2", "CDQ008", "CDQ010", "SDDSRVYR")
  responses <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")
  all_vars <- c(pred, responses)
  
  for (i in seq_along(label_order)) {
    current_label <- label_order[i]
    
    # Define predictors: original predictors + previously predicted labels
    if (i > 1) {
      # Get previous predicted labels
      previous_labels <- label_order[seq_len(i - 1)]
      # Exclude responses (outcomes) from the predictors, add previous labels
      predictors <- paste(c(all_vars[!all_vars %in% c(Y, previous_labels)], previous_labels), collapse = " + ")
      
    } else {
      # For the first model, only use the original predictors (excluding responses)
      predictors <- paste(all_vars[!all_vars %in% Y], collapse = " + ")
    }
    
    # Formula for the current label
    formula <- as.formula(paste(current_label, "~", predictors))
    
    # Train a survey-weighted logistic regression model
    models[[current_label]] <- svyglm(formula, design = survey_design, family = quasibinomial())
  }
  
  return(models)
}

predict_survey_pcc <- function(models, data, order) {
  predictions <- list()  # Initialize list to store predictions
  
  # Loop over each model in the chain
  for (i in seq_along(models)) {
    current_model <- models[[i]]  # Get the current model
    outcome_name <- order[i]      # Get the corresponding outcome name
    
    # Define predictors as all the previous outcomes and the main predictors
    predictor_cols <- names(data)[!names(data) %in% order[i:8]]
    model_data <- data[, predictor_cols]
    
    # Predict for the current model
    predictions[[outcome_name]] <- predict(current_model, newdata = model_data, type = "response")
    
    # Save predictions as a new column in the data
    data[[outcome_name]] <- as.factor(ifelse(predictions[[outcome_name]] > 0.5, 1, 0))
  }
  
  return(predictions)
}





# Generate random permutations of the labels
generate_label_orders <- function(labels, num_orders) {
  all_orders <- replicate(num_orders, sample(labels), simplify = FALSE)
  return(all_orders)
}

# Define outcomes and predictors
outcomes <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")
predictors <- c("DEPR_TOT", "AGE_BIN", "RIAGENDR", "RIDRETH1", "SMQ020", "SMQ040", "HIQ011", 
                "BPQ020", "BMXBMI", "BPQ080", "ALQ130", "MEDDEP", "DMDBORNT", "SDDSRVYR",
                "PAQMV", "DIDTOT", "DUQTOT", "INC3", "DMDEDUC2", "CDQ008", "CDQ010")

# Subset the survey design to include only predictors and outcomes
all_vars_CAD1 <- survey_design_CAD1$variables[, c(predictors, outcomes)]
all_vars_CAD0 <- survey_design_CAD0$variables[, c(predictors, outcomes)]
# Generate label orders
label_orders <- generate_label_orders(outcomes, num_orders = 8)

# Train both models
model_cad1 <- lapply(label_orders, function(order) train_survey_pcc_str(survey_design_CAD1, outcomes, order))
model_cad0 <- lapply(label_orders, function(order) train_survey_pcc_str(survey_design_CAD0, outcomes, order))


# Function to compute variable importance using permutation
compute_variable_importance <- function(survey_design, data, outcomes, predictors, label_orders, ensemble_models, metric_fn, num_permutations = 10) {
  importance_results <- list()
  
  # Get baseline performance before permutation
  baseline_predictions <- list()
  
  for (chain_index in seq_along(ensemble_models)) {
    chain_models <- ensemble_models[[chain_index]]
    chain_order <- label_orders[[chain_index]]
    
    # Get predictions without permutation
    baseline_predictions[[chain_index]] <- predict_survey_pcc(chain_models, data, chain_order)
  }
  
  flattened_predictions <- lapply(seq_along(outcomes), function(idx) {
    outcome <- outcomes[idx]  # Current outcome
    
    # Extract predictions for this outcome across all chains
    chain_predictions <- lapply(baseline_predictions, function(chain) chain[[outcome]])
    
    # Aggregate predictions (mean across chains)
    Reduce("+", chain_predictions) / length(chain_predictions)
  })
  
  # Convert the list to a data frame or matrix for easier handling
  baseline_final_predictions <- do.call(cbind, flattened_predictions)
  
  baseline_metric <- metric_fn(baseline_final_predictions, survey_design)
  
  # Loop through each predictor
  for (predictor in predictors) {
    permuted_scores <- c()
    
    for (perm in seq_len(num_permutations)) {
      permuted_vars <- data
      permuted_vars[[predictor]] <- sample(data[[predictor]])
      print(paste("Permutation #", perm))
      permuted_predictions <- list()
      
      for (chain_index in seq_along(ensemble_models)) {
        chain_models <- ensemble_models[[chain_index]]
        chain_order <- label_orders[[chain_index]]
        
        # Predict using permuted dataset
        permuted_predictions[[chain_index]] <- predict_survey_pcc(chain_models, permuted_vars, chain_order)
      }
      
      
      flattened_predictions <- lapply(seq_along(outcomes), function(idx) {
        outcome <- outcomes[idx]  # Current outcome
        
        # Extract predictions for this outcome across all chains
        chain_predictions <- lapply(permuted_predictions, function(chain) chain[[outcome]])
        
        # Aggregate predictions (mean across chains)
        Reduce("+", chain_predictions) / length(chain_predictions)
      })
      
      # Convert the list to a data frame or matrix for easier handling
      permuted_final_predictions <- do.call(cbind, flattened_predictions)
      
      permuted_metric <- metric_fn(permuted_final_predictions, survey_design)
      
      # Store the drop in performance
      permuted_scores <- c(permuted_scores, permuted_metric - baseline_metric)
    }
    
    # Store mean importance score for the predictor
    importance_results[[predictor]] <- mean(permuted_scores)
  }
  
  return(importance_results)
}

metric_fn <- function(predictions, survey_design) {
  outcomes <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")
  
  # Ensure outcomes are factors and convert them to numeric levels
  true_labels <- survey_design$variables[outcomes]
  true_labels <- as.data.frame(lapply(survey_design$variables[outcomes], function(x) as.numeric(as.character(x))))
  
  # Calculate the Brier score: mean squared error between predicted probabilities and actual class labels
  brier_score <- rowMeans((predictions - true_labels)^2)
  
  # Use svymean to compute the mean Brier score for the weighted survey data
  return(coef(svymean(~ brier_score, survey_design)))
}


# Compute Variable Importance
variable_importance_CAD1 <- compute_variable_importance(survey_design_CAD1, data=all_vars_CAD1, outcomes, predictors, label_orders, model_cad1, metric_fn, num_permutations = 10)
variable_importance_CAD0 <- compute_variable_importance(survey_design_CAD0, data=all_vars_CAD0, outcomes, predictors, label_orders, model_cad0, metric_fn, num_permutations = 10)
# Print variable importance scores
vimp_CAD1 <- unlist(variable_importance_CAD1)
vimp_CAD0 <- unlist(variable_importance_CAD0)



# # Alt choice for metric fn (log loss)
# metric_fn <- function(predictions, survey_design) {
#   # Calculate log-loss (cross-entropy) between predicted and actual values
#   actual_classes <- as.matrix(survey_design$variables[outcomes])
#   epsilon <- 1e-15  # To prevent log(0)
#   predictions <- pmin(pmax(predictions, epsilon), 1 - epsilon)  # Clip predictions between epsilon and 1-epsilon
#   
#   log_loss <- -mean(actual_classes * log(predictions) + (1 - actual_classes) * log(1 - predictions))
#   return(log_loss)
# }




# Conditional VIMP

cond_compute_variable_importance <- function(survey_design, data, outcomes, predictors, label_orders, ensemble_models, metric_fn, num_permutations = 10) {
  importance_results <- matrix(0, nrow = length(predictors), ncol = length(outcomes), dimnames = list(predictors, outcomes))
  
  # Get baseline performance for each outcome separately
  baseline_metrics <- sapply(outcomes, function(outcome) {
    baseline_predictions <- list()
    
    for (chain_index in seq_along(ensemble_models)) {
      chain_models <- ensemble_models[[chain_index]]
      chain_order <- label_orders[[chain_index]]
      baseline_predictions[[chain_index]] <- predict_survey_pcc(chain_models, data, chain_order)
    }
    
    # Compute the final probability for the given outcome
    baseline_final_prediction <- rowMeans(sapply(baseline_predictions, function(chain) chain[[outcome]]))
    
    # Compute baseline performance metric for this outcome
    metric_fn2(baseline_final_prediction, survey_design, outcome)
  })
  
  names(baseline_metrics) <- outcomes
  print(baseline_metrics)
  str(baseline_metrics)
  
  # Loop through each predictor
  for (predictor in predictors) {
    print(paste("Starting predictor:", predictor))
    
    for (outcome in outcomes) {
      permuted_scores <- c()
      
      for (perm in seq_len(num_permutations)) {
        permuted_data <- data
        permuted_data[[predictor]] <- sample(data[[predictor]])
        
        permuted_predictions <- list()
        
        for (chain_index in seq_along(ensemble_models)) {
          chain_models <- ensemble_models[[chain_index]]
          chain_order <- label_orders[[chain_index]]
          permuted_predictions[[chain_index]] <- predict_survey_pcc(chain_models, permuted_data, chain_order)
        }
        
        # Compute final probability for the given outcome
        permuted_final_prediction <- rowMeans(sapply(permuted_predictions, function(chain) chain[[outcome]]))
        
        # Compute performance with permuted predictor
        permuted_metric <- metric_fn2(permuted_final_prediction, survey_design, outcome)
        
        # Store the drop in performance
        permuted_scores <- c(permuted_scores, permuted_metric - baseline_metrics[outcome])
        
        print(paste("Permuted Score:", permuted_scores))
      }
      
      # Store mean importance score for the predictor on this specific outcome
      importance_results[predictor, outcome] <- mean(permuted_scores)
    }
  }
  
  return(importance_results)
}

metric_fn2 <- function(predictions, survey_design, outcome) {
  # Extract observed binary outcomes
  actual_values <- survey_design$variables[[outcome]]
  actual_values <- as.character(actual_values)
  actual_values <- as.numeric(actual_values)
  
  predictions <- as.vector(predictions)
  
  str(predictions)
  str(actual_values)
  
  # Compute squared error
  squared_errors <- (predictions - actual_values)^2
  
  # Compute weighted mean of squared errors using survey weights
  weighted_brier_score <- coef(svymean(~ squared_errors, survey_design))
  
  return(weighted_brier_score)
}




# Compute Conditional VIMP 
cond_variable_importance_CAD1 <- cond_compute_variable_importance(survey_design_CAD1, data=all_vars_CAD1, outcomes, predictors, label_orders, model_cad1, metric_fn2, num_permutations = 10)
cond_variable_importance_CAD0 <- cond_compute_variable_importance(survey_design_CAD0, data=all_vars_CAD0, outcomes, predictors, label_orders, model_cad0, metric_fn2, num_permutations = 10)



library(ggplot2)
library(reshape2)

# Convert matrix to long format
long_data <- melt(vi_dep_cad)

# Plot heatmap
ggplot(long_data, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +  # Customize color scale
  labs(x = "Outcomes", y = "Predictors", fill = "Importance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# VI DEPR_TOT with and without CAD
vi_dep_cad <- rbind(cond_variable_importance_CAD0[1,], cond_variable_importance_CAD1[1,])
rownames(vi_dep_cad) <- c("DEPR_TOT:NOCAD", "DEPR_TOT:CAD")
