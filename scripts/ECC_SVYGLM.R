library(randomForest)
library(randomForestSRC)
library(mice)
library(survey)
library(pROC)
library(splines)

pred <- c("ns(DEPR_TOT,df=3)","AGE_BIN", "RIAGENDR", "RIDRETH1", "SMQ020",
          "SMQ040", "HIQ011", "BPQ020", "ns(BMXBMI,df=3)", "BPQ080",
          "ns(ALQ130,df=3)", "MEDDEP", "DMDBORNT", "PAQMV", "CADTOT", "DIDTOT",
          "DUQTOT", "INC3", "DMDEDUC2", "CDQ008", "CDQ010", "SDDSRVYR")
responses <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")

predictors <- c("DEPR_TOT","AGE_BIN", "RIAGENDR", "RIDRETH1", "SMQ020",
                "SMQ040", "HIQ011", "BPQ020", "BMXBMI", "BPQ080",
                "ALQ130", "MEDDEP", "DMDBORNT", "PAQMV", "CADTOT", "DIDTOT",
                "DUQTOT", "INC3", "DMDEDUC2", "CDQ008", "CDQ010")
wd_subset <- wd[, c(responses, predictors,"SDDSRVYR", "SDMVPSU", "SDMVSTRA", "MEC15YR", "SEQN")]


# Initialize the predictor matrix for MICE
predictorMatrix <- mice::make.predictorMatrix(wd_subset)
predictorMatrix[, !colnames(predictorMatrix) %in% predictors] <- 0
predictorMatrix[!rownames(predictorMatrix) %in% predictors, ] <- 0

# Run MICE
imputations <- mice(wd_subset, seed = 123, maxit = 1, m = 1, predictorMatrix = predictorMatrix, printFlag = TRUE)

# Survey design setup
design <- svydesign(
  id = ~SDMVPSU,
  weights = ~MEC15YR,
  strata = ~SDMVSTRA,
  data = complete(imputations, 1),
  nest = TRUE
)

survey_design <- subset(design, !apply(
  design$variables[, responses], 1,
  function(x) all(x == 0 | is.na(x))
))




train_survey_pcc <- function(survey_design, Y, label_order) {
  models <- list()
  pred <- c("ns(DEPR_TOT,df=3)","AGE_BIN", "RIAGENDR", "RIDRETH1", "SMQ020",
            "SMQ040", "HIQ011", "BPQ020", "ns(BMXBMI,df=3)", "BPQ080",
            "ns(ALQ130,df=3)", "MEDDEP", "DMDBORNT", "PAQMV", "CADTOT", "DIDTOT",
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
      predictors <- paste(c(colnames(survey_design$variables)[!colnames(survey_design$variables) %in% c(Y, previous_labels)], previous_labels), collapse = " + ")
    } else {
      # For the first model, only use the original predictors (excluding responses)
      predictors <- paste(colnames(survey_design$variables)[!colnames(survey_design$variables) %in% Y], collapse = " + ")
    }
    
    # Formula for the current label
    formula <- as.formula(paste(current_label, "~", predictors))
    
    # Train a survey-weighted logistic regression model
    models[[current_label]] <- svyglm(formula, design = survey_design, family = quasibinomial())
  }
  
  return(models)
}

# Function to predict using an ensemble of models with survey data
predict_survey_pcc <- function(models, data, order) {
  predictions <- list()  # Initialize list to store predictions
  
  # Loop over each model in the chain
  for (i in seq_along(models)) {
    current_model <- models[[i]]  # Get the current model
    outcome_name <- order[i]      # Get the corresponding outcome name
    
    # Ensure the model_data includes previous predictions and the current outcome
    model_data <- data$variables[, 
                                          c(setdiff(names(data), order[seq_len(i)]), outcome_name), 
                                          drop = FALSE
    ]
    
    # Debugging: Check that previous predictions exist
    print(paste("Modeling outcome:", outcome_name))
    print("Available predictors:")
    print(names(model_data))
    
    # Predict for the current model
    predictions[[outcome_name]] <- predict(current_model, newdata = model_data, type = "response")
    
    # Debugging: Check predictions before adding to survey design
    print(paste("Predictions for", outcome_name, ":"))
    print(head(predictions[[outcome_name]]))
    
    # Save predictions as a new column in survey_design
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
                "BPQ020", "BMXBMI", "BPQ080", "ALQ130", "MEDDEP", "DMDBORNT", 
                "PAQMV", "CADTOT", "DIDTOT", "DUQTOT", "INC3", "DMDEDUC2", "CDQ008", "CDQ010")

# Subset the survey design to include only predictors and outcomes
all_vars <- survey_design$variables[, c(predictors, outcomes)]

# Generate label orders
label_orders <- generate_label_orders(outcomes, num_orders = 3)

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
print(str(ensemble_predictions))

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

final_binary <- ifelse(final_predictions > 0.5, 1, 0)

# Extract true labels from survey design
true_labels <- survey_design$variables[outcomes]

# Calculate Weighted Hamming Loss
weighted_hamming_loss <- svymean(~ rowMeans(final_binary != true_labels), survey_design)

# Calculate Weighted Subset Accuracy
weighted_subset_accuracy <- svymean(~ I(rowSums(final_binary == true_labels) == ncol(true_labels)), survey_design)

# Print evaluation metrics
print(paste("Weighted Hamming Loss:", weighted_hamming_loss))
print(paste("Weighted Subset Accuracy:", weighted_subset_accuracy))
