library(randomForest)
library(randomForestSRC)
library(mice)
library(survey)
library(pROC)
library(missForest)
library(car)


responses <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")

predictors <- c("DEPR_BIN","AGE_BIN", "RIAGENDR", "RIDRETH1", "SMQ020",
                "SMQ040", "HIQ011", "BPQ020", "BMXBMI", "BPQ080",
                "ALQ130", "MEDDEP", "DMDBORNT", "PAQMV", "CADTOT", "DIDTOT",
                "DUQTOT", "INC3", "DMDEDUC2", "CDQ008", "CDQ010")
wd_subset <- wd[, c(responses, predictors,"SDDSRVYR", "SDMVPSU", "SDMVSTRA", "MEC15YR", "SEQN")]


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
# Perform missing data imputation using missForest
subset <- as.matrix.data.frame(wd_subset)
imputed_data <- missForest(xmis=survey_design$variables, verbose = TRUE)

# Extract the completed dataset
wd_subset_imputed <- imputed_data$ximp  # This is the dataset with imputed values

# Check imputation error (if needed)
print(imputed_data$OOBerror)  # Out-of-bag (OOB) error for numeric/categorical variables

survey_design$variables <- wd_subset_imputed




train_survey_pcc <- function(survey_design, Y, label_order) {
  models <- list()
  vif <- list()
  pred <- c("DEPR_BIN","AGE_BIN", "RIAGENDR", "RIDRETH1", "SMQ020",
            "SMQ040", "HIQ011", "BPQ020", "BMXBMI", "BPQ080",
            "ALQ130", "MEDDEP", "PAQMV", "CADTOT", "DIDTOT",
            "INC3", "DMDEDUC2", "CDQ008", "CDQ010", "SDDSRVYR")
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
    # Function to compute and print high VIF values
    check_vif <- function(model) {
      vif_values <- vif(model)  # Compute VIF
      
      # Extract the normalized VIF (3rd column: GVIF^(1/(2*Df)))
      normalized_vif <- vif_values[, 3]
      
      # Identify variables with high VIF
      high_vif <- normalized_vif[normalized_vif > 5]
      
      if (length(high_vif) > 0) {
        cat("Variables with high VIF (>5):\n")
        
        # Create a dataframe for better readability
        high_vif_df <- data.frame(
          Variable = names(high_vif),
          VIF_Value = high_vif
        )
        
        print(high_vif_df, row.names = FALSE)  # Print without row indices
      } else {
        cat("No multicollinearity issues detected (all VIF < 5).\n")
      }
      
    }
    
    
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
predictors <- c("DEPR_BIN", "AGE_BIN", "RIAGENDR", "RIDRETH1", "SMQ020", "SMQ040", "HIQ011", 
                "BPQ020", "BMXBMI", "BPQ080", "ALQ130", "MEDDEP", "SDDSRVYR",
                "PAQMV", "CADTOT", "DIDTOT", "INC3", "DMDEDUC2", "CDQ008", "CDQ010")

# Subset the survey design to include only predictors and outcomes
all_vars <- survey_design$variables[, c(predictors, outcomes)]

# Generate label orders
label_orders <- generate_label_orders(outcomes, num_orders = 8)

# Train PCC ensemble
ensemble_models <- lapply(label_orders, function(order) train_survey_pcc(survey_design, outcomes, order))

# Predict with the ensemble
# Loop through each chain and make predictions for each arrangement
ensemble_predictions <- list()  # List to store the predictions for all chains

# Iterate over each chain arrangement
for (chain_index in seq_along(ensemble_models)) {
  # Get the models for the current chain and the order of outcomes
  chain_models <- ensemble_models[[chain_index]]
  chain_order <- label_orders[[chain_index]]
  
  # Debugging: Check the chain index and order
  print(paste("Processing Chain", chain_index, "with order:", paste(chain_order, collapse = ", ")))
  
  # Get the predictions for the current chain
  ensemble_predictions[[chain_index]] <- predict_survey_pcc(chain_models, all_vars, chain_order)
}

# Debugging: Check the final structure of ensemble_predictions
print("Ensemble predictions completed:")

# Combine predictions for each outcome across chains
flattened_predictions <- lapply(seq_along(outcomes), function(idx) {
  outcome <- outcomes[idx]  # Current outcome
  
  # Extract predictions for this outcome across all chains
  chain_predictions <- lapply(ensemble_predictions, function(chain) chain[[outcome]])
  
  # Aggregate predictions (mean across chains)
  Reduce("+", chain_predictions) / length(chain_predictions)
})

# Convert the list to a data frame or matrix for easier handling
final_predictions <- do.call(cbind, flattened_predictions)
colnames(final_predictions) <- outcomes

final_binary <- ifelse(final_predictions > 0.33, 1, 0)

# Extract true labels from survey design
true_labels <- survey_design$variables[outcomes]
true_labels <- as.data.frame(lapply(survey_design$variables[outcomes], function(x) as.numeric(as.character(x))))


# Compute weighted TP, FP, and FN for each outcome
weighted_tp <- svymean(~ rowMeans((final_binary == 1) & (true_labels == 1)), survey_design)
weighted_fp <- svymean(~ rowMeans((final_binary == 1) & (true_labels == 0)), survey_design)
weighted_fn <- svymean(~ rowMeans((final_binary == 0) & (true_labels == 1)), survey_design)

# Extract coefficient values
tp_value <- coef(weighted_tp)
fp_value <- coef(weighted_fp)
fn_value <- coef(weighted_fn)

# Compute precision, recall, and F1-score
precision <- tp_value / (tp_value + fp_value)
recall <- tp_value / (tp_value + fn_value)
f1_score <- 2 * (precision * recall) / (precision + recall)

# Print results
print(paste("Survey-weighted Precision:", precision))
print(paste("Survey-weighted Recall:", recall))
print(paste("Survey-weighted F1-score:", f1_score))



# Print evaluation metrics
weighted_brier_score <- svymean(~ rowMeans((final_predictions - true_labels)^2), survey_design)
weighted_prob_hamming_loss <- svymean(~ rowMeans(abs(final_predictions - true_labels)), survey_design)
expected_subset_accuracy <- svymean(~ rowMeans(final_binary == true_labels), survey_design)
expected_subset_accuracy_value <- coef(expected_subset_accuracy)
print(paste("Expected Subset Accuracy:", expected_subset_accuracy_value))
print(paste("Probability-Weighted Hamming Loss:", weighted_prob_hamming_loss))
print(paste("Weighted Brier Score:", weighted_brier_score))



# Function to compute variable importance using permutation
compute_variable_importance <- function(survey_design, data, outcomes, predictors, label_orders, ensemble_models, metric_fn, num_permutations = 50) {
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
variable_importance <- compute_variable_importance(survey_design, data=all_vars, outcomes, predictors, label_orders, ensemble_models, metric_fn, num_permutations = 10)

# Print variable importance scores
vimp <- unlist(variable_importance)
print(vimp)


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
cond_variable_importance <- cond_compute_variable_importance(survey_design, data=all_vars, outcomes, predictors, label_orders, ensemble_models, metric_fn2, num_permutations = 10)
cond_vimp <- unlist(cond_variable_importance)
print(cond_vimp)



library(ggplot2)
library(reshape2)

# Convert matrix to long format
long_data <- melt(cond_variable_importance)

# Plot heatmap
ggplot(long_data, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +  # Customize color scale
  labs(x = "Outcomes", y = "Predictors", fill = "Importance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

