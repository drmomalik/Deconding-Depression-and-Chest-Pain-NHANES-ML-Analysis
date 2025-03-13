# Assuming the following:
# chain_models: Your ensemble classification chain model (list of models)
# data: The original dataset with all variables
# chain_order: The order of the responses in the chain
# DEPR_TOT is the predictor you want to vary

outcomes <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")
predictors <- c("DEPR_TOT", "AGE_BIN", "RIAGENDR", "RIDRETH1", "SMQ020", "SMQ040", "HIQ011", 
                "BPQ020", "BMXBMI", "BPQ080", "ALQ130", "MEDDEP", "DMDBORNT", "SDDSRVYR",
                "PAQMV", "CADTOT", "DIDTOT", "DUQTOT", "INC3", "DMDEDUC2", "CDQ008", "CDQ010")
survey_design$variables <- wd_subset_imputed
data <- survey_design$variables[, c(predictors, outcomes)]

# Function to calculate the most common level for categorical variables
get_mode <- function(x) {
  uniq_x <- unique(x)
  uniq_x[which.max(tabulate(match(x, uniq_x)))]  # Find the most frequent value
}

# 1. Initialize an empty list to store the calculated values
medians <- vector("list", length = ncol(data))

# 2. Loop over each column and calculate median or mode, preserving data types
for (i in seq_along(data)) {
  if (is.numeric(data[[i]])) {
    # For numeric variables, calculate the median
    medians[[i]] <- median(data[[i]], na.rm = TRUE)
  } else if (is.factor(data[[i]])) {
    # For factor variables, calculate the mode
    medians[[i]] <- get_mode(data[[i]])
  } else if (is.character(data[[i]])) {
    # For character variables, calculate the mode (if applicable)
    medians[[i]] <- get_mode(data[[i]])
  }
}

# 3. Convert the list of medians back to a dataframe while preserving data types
medians_df <- as.data.frame(medians)

# 4. Set the variable names to match the original column names
colnames(medians_df) <- colnames(data)

# 5. Set DEPR_TOT to NA since we want it to vary
medians_df$DEPR_TOT <- NA  # Leave DEPR_TOT to vary
medians_df$DEPR_TOT <- as.numeric(medians_df$DEPR_TOT)

# View the result
head(medians_df)


# 2. Generate a sequence of DEPR_TOT values (min to max)
depr_tot_seq <- seq(from = floor(min(data$DEPR_TOT, na.rm = TRUE)), 
                    to = ceiling(max(data$DEPR_TOT, na.rm = TRUE)), 
                    by = 1)

# 3. Create a dataframe to store the predictions
predictions_df <- data.frame()

# 4. Loop over the DEPR_TOT sequence and predict for each value
for (depr_value in depr_tot_seq) {
  
  # Set the current value of DEPR_TOT
  medians_df["DEPR_TOT"] <- depr_value
  
  # Create a new dataset with the current DEPR_TOT and constant other variables
  temp_data <- medians_df
  
  # Predict using the `predict_survey_pcc` function
  for (chain_index in seq_along(ensemble_models)) {
    # Get the models for the current chain and the order of outcomes
    chain_models <- ensemble_models[[chain_index]]
    chain_order <- label_orders[[chain_index]]
    
    # Debugging: Check the chain index and order
    print(paste("Processing Chain", chain_index, "with order:", paste(chain_order, collapse = ", ")))
    
    # Get the predictions for the current chain
    ensemble_predictions[[chain_index]] <- predict_survey_pcc(chain_models, temp_data, chain_order)
  }
  
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
  print(final_predictions)
  preds <- final_predictions
  
  # Store the predictions (probabilities) for each outcome
  predictions_df <- rbind(predictions_df, preds)
}

rownames(predictions_df) <- depr_tot_seq

# Ensure DEPR_TOT is included as a column, rather than row names
predictions_df$DEPR_TOT <- rownames(predictions_df)

# Reshape the data from wide to long format
library(tidyr)
long_df <- predictions_df %>%
  pivot_longer(cols = -DEPR_TOT, 
               names_to = "Outcome", 
               values_to = "Probability")

# Convert DEPR_TOT to numeric (if necessary)
long_df$DEPR_TOT <- as.numeric(long_df$DEPR_TOT)

# Plot the probabilities of each outcome as DEPR_TOT increases
library(ggplot2)
ggplot(long_df, aes(x = DEPR_TOT, y = Probability, color = Outcome)) +
  geom_line() + 
  labs(title = "Probability of Outcomes as DEPR_TOT Increases",
       x = "DEPR_TOT",
       y = "Probability") +
  theme_minimal() +
  theme(legend.title = element_blank())