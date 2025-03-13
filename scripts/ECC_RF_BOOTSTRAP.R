library(randomForest)
library(randomForestSRC)
library(ranger)
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

rep_design <- as.svrepdesign(survey_design, type = "bootstrap", replicates = 10)

# Stabilize weights
repweights <- weights(rep_design, type = "analysis") / mean(weights(rep_design, type = "analysis"))
pweights <- weights(rep_design, type = "sampling") / mean(weights(rep_design, type = "sampling"))

data <- rep_design$variables

train_survey_pcc <- function(data, outcomes, label_order) {
  ecc <- list()
  for (o in seq_along(label_order)) {
    print(paste("On chain #", o))
    order <- label_order[[o]]
    pred <- c("DEPR_BIN","AGE_BIN", "RIAGENDR", "RIDRETH1", "SMQ020",
              "SMQ040", "HIQ011", "BPQ020", "BMXBMI", "BPQ080",
              "ALQ130", "MEDDEP", "PAQMV", "CADTOT", "DIDTOT","DMDBORNT",
              "INC3", "DMDEDUC2", "CDQ008", "CDQ010", "SDDSRVYR", "DUQTOT")
    responses <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")
    all_vars <- c(pred, responses)
    
    models <- list()
    
    for (i in seq_along(order)) {
      current_label <- order[i]
      print(paste("Link #", i, "of", length(order)))
      
      # Define predictors: original predictors + previously predicted labels
      if (i > 1) {
        # Get previous predicted labels
        previous_labels <- order[seq_len(i - 1)]
        # Exclude responses (outcomes) from the predictors, add previous labels
        predictors <- paste(c(all_vars[!all_vars %in% c(outcomes, previous_labels)], previous_labels), collapse = " + ")
        
      } else {
        # For the first model, only use the original predictors (excluding responses)
        predictors <- paste(all_vars[!all_vars %in% outcomes], collapse = " + ")
      }
      
      # Formula for the current label
      formula <- as.formula(paste(current_label, "~", predictors))
      
      # Train a survey-weighted logistic regression model
      models[[current_label]] <- imbalanced(formula, data=data, case.weights = data$weights, 
                                        method="rfq", splitrule = "auc", num.trees = 1000, samptype="swr",
                                        nodesize=1, fast=TRUE, forest=TRUE)
     
      }
    
    ecc[[o]] <- models
    }
    
    return(ecc)
  }
  
  predict_survey_pcc <- function(models, data, order) {
    predictions <- list()  # Initialize list to store predictions
    
    # Loop over each model in the chain
    for (i in seq_along(models)) {
      current_model <- models[[i]]  # Get the current model
      outcome_name <- order[i]      # Get the corresponding outcome name
      
      # Define predictors as all the previous outcomes and the main predictors
      if (i==8) {
        model_data <- data
      }
      else {
        predictor_cols <- names(data)[!names(data) %in% order[(i+1):8]]
        model_data <- data[, predictor_cols]
      }
     
      
      # Predict for the current model
      pred_object <- predict(current_model, newdata = model_data, outcome="test")
      predictions[[outcome_name]] <- pred_object$class.oob

      # Save predictions as a new column in the data
      data[[outcome_name]] <- predictions[[outcome_name]]
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
                "BPQ020", "BMXBMI", "BPQ080", "ALQ130", "MEDDEP", "SDDSRVYR",
                "PAQMV", "DUQTOT", "DMDBORNT", "CADTOT", "DIDTOT", "INC3", "DMDEDUC2", "CDQ008", "CDQ010")



# Generate label orders
label_orders <- generate_label_orders(outcomes, num_orders = 5)

# Train PCC ensemble
ensemble_models <- list()
for (w in 1:ncol(train_repweights)) {
  print(paste("On bootstrap sample #:", w))
  train_data$weights <- train_repweights[,w]
  ensemble_models[[w]] <- train_survey_pcc(train_data, outcomes, label_orders)
}


replicate_results <- data.frame(rep = integer(),
                                expected_subset_accuracy = numeric(),
                                f1 = numeric(),
                                brier_score = numeric(),
                                prob_hamming_loss = numeric(),
                                TNR = numeric(),
                                TPR = numeric(),
                                G_Mean = numeric(),
                                stringsAsFactors = FALSE)

for(ecc in seq_along(ensemble_models)) {
  
  model <- ensemble_models[[ecc]]
  # Predict with the ensemble
  # Loop through each chain and make predictions for each arrangement
  ensemble_predictions <- list()  # List to store the predictions for all chains
  
  # Iterate over each chain arrangement
  for (chain_index in seq_along(model)) {
    # Get the models for the current chain and the order of outcomes
    chain_models <- model[[chain_index]]
    chain_order <- label_orders[[chain_index]]
    
    # Debugging: Check the chain index and order
    print(paste("Processing Chain", chain_index, "with order:", paste(chain_order, collapse = ", ")))
    
    # Get the predictions for the current chain
    ensemble_predictions[[chain_index]] <- predict_survey_pcc(chain_models, test_data, chain_order)
  }
  
  # Debugging: Check the final structure of ensemble_predictions
  print("Ensemble predictions completed:")

  # Combine predictions for each outcome across chains
  flattened_predictions <- lapply(seq_along(outcomes), function(idx) {
    outcome <- outcomes[idx]  # Current outcome
    
    # Extract predictions for this outcome across all chains
    chain_predictions <- lapply(ensemble_predictions, function(chain) chain[[outcome]])
    # Convert the factors in each vector to numeric (0 or 1)
    numeric_predictions <- lapply(chain_predictions, function(x) as.numeric(as.character(x)))
    
    # Calculate the mean prediction across all the chains
    Reduce("+", numeric_predictions) / length(numeric_predictions)
  })

  # Convert the list to a data frame or matrix for easier handling
  final_predictions <- do.call(cbind, flattened_predictions)
  colnames(final_predictions) <- outcomes

  # Extract true labels from survey design
  true_labels <- test_rep$variables[outcomes]
  true_labels <-as.matrix(sapply(true_labels, function(x) as.numeric(as.character(x))))
  
  # Compute prevalence π (P(Y = 1)) for each outcome in true labels
  prevalence <- colMeans(true_labels)  # Proportion of minority class (1s)
  
  # Apply quantile classifier threshold for each outcome
  final_binary_qstar <- final_predictions
  for (i in seq_along(outcomes)) {
    final_binary_qstar[, i] <- ifelse(final_predictions[, i] >= prevalence[i], 1, 0)
  }
  
  # Initialize performance storage
  TNR_values <- numeric(length(outcomes))
  TPR_values <- numeric(length(outcomes))
  G_Mean_values <- numeric(length(outcomes))
  
  # Compute TNR, TPR, and G-Mean for each outcome
  for (i in seq_along(outcomes)) {
    cm <- table(Predicted = final_binary_qstar[, i], Actual = true_labels[, i])
    
    # Extract values from confusion matrix
    tn <- ifelse("0" %in% rownames(cm) & "0" %in% colnames(cm), cm["0", "0"], 0)
    fp <- ifelse("1" %in% rownames(cm) & "0" %in% colnames(cm), cm["1", "0"], 0)
    fn <- ifelse("0" %in% rownames(cm) & "1" %in% colnames(cm), cm["0", "1"], 0)
    tp <- ifelse("1" %in% rownames(cm) & "1" %in% colnames(cm), cm["1", "1"], 0)
    
    # Compute TNR, TPR
    TNR_values[i] <- tn / (tn + fp)
    TPR_values[i] <- tp / (tp + fn)
    
    # Compute G-Mean
    G_Mean_values[i] <- sqrt(TNR_values[i] * TPR_values[i])
  }
  
  # Store results in a data frame
  performance_results <- data.frame(
    Outcome = outcomes,
    TNR = round(TNR_values, 4),
    TPR = round(TPR_values, 4),
    G_Mean = round(G_Mean_values, 4)
  )
  
  
  # ----- Metric 1: Expected Subset Accuracy -----
  # This version calculates, for each observation, the average proportion of labels
  # that are correctly predicted.
  expected_subset_accuracy <- mean(apply(final_binary_qstar == true_labels, 1, mean))
  
  # ----- Metric 2: Micro-Averaged F1 Score -----
  # Calculate overall TP, FP, and FN across all outcomes.
  TP <- sum(final_binary_qstar == 1 & true_labels == 1)
  FP <- sum(final_binary_qstar == 1 & true_labels == 0)
  FN <- sum(final_binary_qstar == 0 & true_labels == 1)
  
  # Compute precision and recall; include a check to avoid division by zero.
  precision <- ifelse((TP + FP) > 0, TP / (TP + FP), NA)
  recall    <- ifelse((TP + FN) > 0, TP / (TP + FN), NA)
  f1 <- ifelse((precision + recall) > 0, 2 * precision * recall / (precision + recall), NA)
  
  # ----- Metric 3: Brier Score -----
  # The Brier score is the average squared difference between predicted probabilities and true labels.
  brier_score <- mean((final_predictions - true_labels)^2)
  
  # ----- Metric 4: Probabilistic Hamming Loss -----
  # This is the average absolute difference between the predicted probabilities and the true labels.
  prob_hamming_loss <- mean(abs(final_predictions - true_labels))
  
  # Store the metrics for the current replicate in the data frame.
  replicate_results <- rbind(replicate_results,
                             data.frame(rep = ecc,
                                        expected_subset_accuracy = expected_subset_accuracy,
                                        f1 = f1,
                                        brier_score = brier_score,
                                        prob_hamming_loss = prob_hamming_loss,
                                        TNR = mean(round(TNR_values, 4)),
                                        TPR = mean(round(TPR_values, 4)),
                                        G_Mean = mean(round(G_Mean_values, 4))))
}

# After looping over all replicates, you can average the metrics across replicates.
avg_expected_subset_accuracy <- mean(replicate_results$expected_subset_accuracy, na.rm = TRUE)
avg_f1 <- mean(replicate_results$f1, na.rm = TRUE)
avg_brier_score <- mean(replicate_results$brier_score, na.rm = TRUE)
avg_prob_hamming_loss <- mean(replicate_results$prob_hamming_loss, na.rm = TRUE)

# Print out the averaged metrics
print(paste("Average Expected Subset Accuracy:", avg_expected_subset_accuracy))
print(paste("Average F1 Score:", avg_f1))
print(paste("Average Brier Score:", avg_brier_score))
print(paste("Average Probabilistic Hamming Loss:", avg_prob_hamming_loss))


  
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

