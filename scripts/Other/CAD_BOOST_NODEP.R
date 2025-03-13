library(survey)
library(missForest)
library(splitTools)
library(rsample)
library(xgboost)
library(dplyr)
library(ggplot2)

# Define outcome variables
responses <- c("CAD")
wd$TCP <- ifelse(wd$CDQ009D ==1 | wd$CDQ009E ==1|wd$CDQ009F ==1, 1, 0)
wd$ATCP <- ifelse(!wd$TCP, 1, 0) 
predictors <- c("AGE_BIN", "RIAGENDR", "RIDRETH1", "SMQ020",
                "SMQ040", "HIQ011", "BPQ020", "BMXBMI", "BPQ080",
                "ALQ130", "PAQMV", "DIDTOT","DMDBORNT",
                "INC3", "DMDEDUC2", "SDDSRVYR", "DUQTOT", "CDQ008",
                "CDQ010",
                regions)
regions <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")
wd_subset <- wd[, c(responses, predictors, "CDQ001", "DEPR_NEW", "MEDDEP", "HUQ090", "SDMVPSU", "SDMVSTRA", "MEC15YR", "SEQN")]


# Extract the data from the svrepdesign object
design <- svydesign(
  id = ~SDMVPSU,
  weights = ~MEC15YR,
  strata = ~SDMVSTRA,
  data = wd_subset,
  nest = TRUE
)

survey_design <- subset(design, design$variables["CDQ001"] == 1 & design$variables["DEPR_NEW"] == 0 & design$variables["MEDDEP"] == 0
                        & design$variables["HUQ090"] == 2)

strata <- survey_design$variables$SDMVSTRA
clusters <- survey_design$variables$SDMVPSU


survey_design$variables$SDMVSTRA <- NULL
survey_design$variables$SDMVPSU <- NULL
survey_design$variables$SEQN <- NULL
survey_design$variables$CDQ001 <- NULL


imputed_data <- missForest(xmis=survey_design$variables, verbose = TRUE)
# Extract the completed dataset
wd_subset_imputed <- imputed_data$ximp  # This is the dataset with imputed values
print(imputed_data$OOBerror)


survey_design$variables <- data <- wd_subset_imputed
survey_design$variables$SDMVSTRA <- strata
survey_design$variables$SDMVPSU <- clusters 

# create full data replicate samples
rep <- as.svrepdesign(survey_design, type="bootstrap", replicates=100)

repweights <- weights(rep, type="analysis")
repweights <- repweights/mean(repweights)

# Full data
data <- survey_design$variables

# Define predictors and label
predictors <- c("AGE_BIN", "RIAGENDR", "RIDRETH1", "SMQ020",
                "SMQ040", "HIQ011", "BPQ020", "BMXBMI", "BPQ080",
                "ALQ130", "PAQMV", "DIDTOT","DMDBORNT",
                "INC3", "DMDEDUC2", "SDDSRVYR", "DUQTOT", "CDQ008",
                "CDQ010",
                regions)
current_label <- "CAD"
formula <- as.formula(paste(current_label, "~", paste(c(predictors), collapse="+")))

data$CAD <- factor(data$CAD)
levels(data$CAD) <- c(0,1)


### Create list of dataframes with each replicate weight attached
all_data <- list()
for (i in 1:ncol(repweights)) {
  df <- data
  df$weights <- repweights[,i]
  all_data[[i]] <- df
}

# In-bag and Out-of-bag data 
ib_data <- list()
samp_list <- list()
df_sampled <- list()
feat_list <- list()
metrics_kfold <- data.frame(Model_Iteration = integer(0), 
                         IF_AUC = numeric(0), 
                         OOF_AUC = numeric(0))
all_boost0 <- list()
for (i in 1:length(all_data)) {
  print(paste("On rep:", i))
  df <- all_data[[i]]
  
  # Sample size (downsample to match smaller depression cohort)
  n <- 1118
  
  # Sample n rows from the dataframe df, stratified by the variable CAD
  df_sampled <- df %>%
    group_by(CAD) %>%
    sample_frac(size = n / nrow(df), replace = TRUE) %>%
    ungroup()
  
  samp_list[[i]] <- df_sampled
  ib_data[[i]] <- df_sampled[df_sampled$weights>0,]
  ob_data[[i]] <- df_sampled[df_sampled$weights==0,]
    
  
  # Identify factor and numeric predictors
  factor_predictors <- predictors[sapply(df_sampled[, predictors], is.factor)]
  numeric_predictors <- predictors[!sapply(df_sampled[, predictors], is.factor)]
  
  # One-hot encode factor predictors
  if (length(factor_predictors) > 0) {
    factor_features <- model.matrix(~ . - 1, data = df_sampled[, factor_predictors])
  } else {
    factor_features <- NULL
  }
  
  # Keep numeric predictors unchanged
  numeric_features <- as.matrix(df_sampled[, numeric_predictors])
  
  # Combine both
  features <- if (!is.null(factor_features)) {
    cbind(numeric_features, factor_features)
  } else {
    numeric_features  # Only numeric predictors case
  }
  
  feat_list[[i]] <- features
  
  # Convert target variable to numeric
  labels <- as.numeric(df_sampled$CAD)-1  # Ensure it's 0/1
  weights <- df_sampled$weights  # Sample weights
  
  # Convert to XGBoost DMatrix
  dtrain <- xgb.DMatrix(data = features, label = labels, weight = weights)
  
  # Set parameters for the XGBoost model
  params <- list(
    objective = "binary:logistic",    # Binary classification
    eval_metric = "auc",           # Default evaluation metric
    max_depth = 6,                    # Example hyperparameter
    eta = 0.1                        # Learning rate
  )
  
  # Perform cross-validation with custom evaluation function
  cv_results <- xgb.cv(
    params = params,
    data = dtrain,
    weights = weights,
    metrics = list("logloss", "auc"),
    nfold = 10,                        # 5-fold cross-validation
    nrounds = 500,                    # Number of boosting rounds
    verbose = 1,                      # Print progress
    print_every_n = 50,                # Print every 10 rounds
    early_stopping_rounds = 50  # Stop after 50 rounds of no improvement
  )

  # Extract best model and store
  best_model <- cv_results$best_iteration
  final_model <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = best_model
  )
  all_boost0[[i]] <- final_model 
  
  # Extract AUC In-fold and OO-fold
  train_auc <- cv_results$evaluation_log$train_auc_mean[[best_model]]
  test_auc <- cv_results$evaluation_log$test_auc_mean[[best_model]]
  
  # Add the current model iteration's results to the dataframe
  metrics_kfold <- rbind(metrics_kfold, data.frame(Model_Iteration = i, 
                                             IF_AUC = train_auc, 
                                             OOF_AUC = test_auc))
}


# Average IF_AUC and OOF_AUC
print(paste("IF_AUC:", mean(metrics_kfold$IF_AUC)))
print(paste("OOF_AUC:", mean(metrics_kfold$OOF_AUC)))


results_df <- data.frame(
  IB_AUC = numeric(), OOB_AUC = numeric(),
  IB_F1 = numeric(), OOB_F1 = numeric(),
  IB_Spec = numeric(), OOB_Spec = numeric(),
  IB_Sens = numeric(), OOB_Sens = numeric(),
  IB_Brier = numeric(), OOB_Brier = numeric()
)

# Calculate other metrics using IB and OOB datasets 
  
for (i in 1:ncol(repweights)) {
  ib_df <- ib_data[[i]]
  ob_df <- ob_data[[i]]
  final_model <- all_boost0[[i]]
    
  #### IB ERROR
  
  # Identify factor and numeric predictors
  factor_predictors <- predictors[sapply(ib_df[, predictors], is.factor)]
  numeric_predictors <- predictors[!sapply(ib_df[, predictors], is.factor)]
  
  # One-hot encode factor predictors
  if (length(factor_predictors) > 0) {
    factor_features <- model.matrix(~ . - 1, data = ib_df[, factor_predictors])
  } else {
    factor_features <- NULL
  }
  
  # Keep numeric predictors unchanged
  numeric_features <- as.matrix(ib_df[, numeric_predictors])
  
  # Combine both
  new_features <- if (!is.null(factor_features)) {
    cbind(numeric_features, factor_features)
  } else {
    numeric_features  # Only numeric predictors case
  }
  
  ib_new_labels <- as.numeric(ib_df$CAD)-1
  
  # Convert new features to xgb.DMatrix
  new_dmatrix <- xgb.DMatrix(data = new_features, label = ib_new_labels)
  
  # Generate predicted probabilities
  ib_predicted_probs <- predict(final_model, new_dmatrix)
  
  # Convert probabilities to predicted classes
  ib_predicted_classes <- ifelse(ib_predicted_probs > 0.5, 1, 0)
  
  # Calculate AUC
  ib_roc_curve <- roc(ib_new_labels, ib_predicted_probs)
  ib_auc_value <- pROC::auc(ib_roc_curve)

  # Calculate F1 Score
  ib_f1_value <- F1_Score(y_true = ib_new_labels, y_pred = ib_predicted_classes)

  # Calculate Specificity and Sensitivity
  conf_matrix <- table(ib_new_labels, ib_predicted_classes)
  ib_sensitivity <- conf_matrix[1, 1] / (conf_matrix[1, 1] + conf_matrix[1, 2])
  ib_specificity <- conf_matrix[2, 2] / (conf_matrix[2, 1] + conf_matrix[2, 2])

  # Calculate Brier Score
  ib_brier_score <- mean((ib_predicted_probs - ib_new_labels)^2)
  
  ### OOB Error
  
  # Identify factor and numeric predictors
  factor_predictors <- predictors[sapply(ob_df[, predictors], is.factor)]
  numeric_predictors <- predictors[!sapply(ob_df[, predictors], is.factor)]
  
  # One-hot encode factor predictors
  if (length(factor_predictors) > 0) {
    factor_features <- model.matrix(~ . - 1, data = ob_df[, factor_predictors])
  } else {
    factor_features <- NULL
  }
  
  # Keep numeric predictors unchanged
  numeric_features <- as.matrix(ob_df[, numeric_predictors])
  
  # Combine both
  new_features <- if (!is.null(factor_features)) {
    cbind(numeric_features, factor_features)
  } else {
    numeric_features  # Only numeric predictors case
  }
  
  oob_new_labels <- as.numeric(ob_df$CAD)-1
  
  # Convert new features to xgb.DMatrix
  new_dmatrix <- xgb.DMatrix(data = new_features, label = oob_new_labels)
  
  # Generate predicted probabilities
  oob_predicted_probs <- predict(final_model, new_dmatrix)
  
  # Convert probabilities to predicted classes
  oob_predicted_classes <- ifelse(oob_predicted_probs > 0.5, 1, 0)
  
  # Calculate AUC
  oob_roc_curve <- roc(oob_new_labels, oob_predicted_probs)
  oob_auc_value <- pROC::auc(oob_roc_curve)
  
  # Calculate F1 Score
  oob_f1_value <- F1_Score(y_true = oob_new_labels, y_pred = oob_predicted_classes)
  
  # Calculate Specificity and Sensitivity
  conf_matrix <- table(oob_new_labels, oob_predicted_classes)
  oob_sensitivity <- conf_matrix[1, 1] / (conf_matrix[1, 1] + conf_matrix[1, 2])
  oob_specificity <- conf_matrix[2, 2] / (conf_matrix[2, 1] + conf_matrix[2, 2])
  
  # Calculate Brier Score
  oob_brier_score <- mean((ib_predicted_probs - ib_new_labels)^2)
  
  
  ## Store results
  results_df <- rbind(results_df, data.frame(
    IB_AUC=ib_auc_value, OOB_AUC=oob_auc_value, 
    IB_F1=ib_f1_value, OOB_F1=oob_f1_value, 
    IB_Spec=ib_specificity, OOB_Spec=oob_specificity, 
    IB_Sens=ib_sensitivity, OOB_Sens=oob_sensitivity,
    IB_Brier=ib_brier_score, OOB_Brier=oob_brier_score
  ))
}  

# Compute the 0.632 Bootstrap Estimates
results_df <- results_df %>%
  mutate(
    `0.632_AUC` = 0.368 * IB_AUC + 0.632 * OOB_AUC,
    `0.632_F1` = 0.368 * IB_F1 + 0.632 * OOB_F1,
    `0.632_Spec` = 0.368 * IB_Spec + 0.632 * OOB_Spec,
    `0.632_Sens` = 0.368 * IB_Sens + 0.632 * OOB_Sens,
    `0.632_Brier` = 0.368 * IB_Brier + 0.632 * OOB_Brier
  )

# Create a final dataframe with only the 0.632 estimates
bootstrap_results <- results_df %>%
  dplyr::select(`0.632_AUC`, `0.632_F1`, `0.632_Spec`, `0.632_Sens`, `0.632_Brier`)

# Calculate summary statistics for each metric
summary_stats <- sapply(bootstrap_results, function(x) {
  c(
    mean = mean(x),
    sd = sd(x),
    `2.5%` = quantile(x, 0.025),
    `25%` = quantile(x, 0.25),
    `50%` = quantile(x, 0.5),
    `75%` = quantile(x, 0.75),
    `97.5%` = quantile(x, 0.975)
  )
})

# Transpose and convert to dataframe
summary_df <- as.data.frame(t(summary_stats))

# Clean metric names (remove "0.632_" prefix)
rownames(summary_df) <- gsub("0.632_", "", rownames(summary_df))

# Format column names
colnames(summary_df) <- c("Mean", "SD", "2.5%", "25%", "50%", "75%", "97.5%")

# Show the resulting dataframe
summary_df



### Extract variable importance 

# Get all feature names from original data
feature_names <- colnames(features)

# Initialize list to store importance data from each model
importance_list <- list()

# Loop through all models in the meta-model list
for (i in seq_along(all_boost0)) {
  # Get current model
  model <- all_boost0[[i]]
  
  # Extract importance matrix
  imp_matrix <- xgboost::xgb.importance(
    model = model,
    feature_names = feature_names
  )
  
  # Create full feature importance dataframe with zeros for missing features
  imp_full <- data.frame(Feature = feature_names) |>
    dplyr::left_join(
      imp_matrix[, c("Feature", "Gain")],  # Use 'Gain' for importance
      by = "Feature"
    ) |>
    dplyr::mutate(
      Importance = ifelse(is.na(.data$Gain), 0, .data$Gain),
      Model = i
    ) |>
    dplyr::select(-Gain)
  
  importance_list[[i]] <- imp_full
}

# Combine all results into single dataframe
combined_importance <- dplyr::bind_rows(importance_list)

# Calculate distribution metrics
distribution_metrics <- combined_importance |>
  dplyr::group_by(Feature) |>
  dplyr::summarise(
    mean = mean(Importance),
    sd = sd(Importance),
    q2.5 = quantile(Importance, 0.025),
    q25 = quantile(Importance, 0.25),
    median = median(Importance),
    q75 = quantile(Importance, 0.75),
    q97.5 = quantile(Importance, 0.975)
  ) |>
  dplyr::arrange(-mean)

# Create box-whisker plot with custom quantiles
ggplot(distribution_metrics, aes(x = reorder(Feature, mean))) +
  geom_boxplot(
    aes(
      ymin = q2.5,
      lower = q25,
      middle = median,
      upper = q75,
      ymax = q97.5
    ),
    stat = "identity",
    fill = "lightblue",
    color = "darkblue"
  ) +
  coord_flip() +
  labs(x = "Feature", y = "Importance Score",
       title = "Feature Importance Distribution Across Models") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))



#### ALE 

bin_vars <- c("CDQ0081", "CDQ0101", "CDQ009A1", "CDQ009B1", "CDQ009C1", "CDQ009D1",
              "CDQ009E1", "CDQ009F1", "CDQ009G1", "CDQ009H1")

ale_val <- list(
  "CDQ008" = vector(length = length(all_boost0)),
  "CDQ010" = vector(length = length(all_boost0)),
  "CDQ009A" = vector(length = length(all_boost0)),
  "CDQ009B" = vector(length = length(all_boost0)),
  "CDQ009C" = vector(length = length(all_boost0)),
  "CDQ009D" = vector(length = length(all_boost0)),
  "CDQ009E" = vector(length = length(all_boost0)),
  "CDQ009F" = vector(length = length(all_boost0)),
  "CDQ009G" = vector(length = length(all_boost0)),
  "CDQ009H" = vector(length = length(all_boost0))
)

for (i in 1:ncol(repweights)) {
  print(paste("On model #:", i))
  final_model <- all_boost0[[i]]
  feat <- as.data.frame(feat_list[[i]])
  
  # Define a prediction function for XGBoost
  predict_function <- function(model, newdata) {
    newdata_matrix <- xgboost::xgb.DMatrix(data.matrix(newdata))
    predict(model, newdata_matrix)
  }
  
  # Create a Predictor object
  predictor <- Predictor$new(
    model = final_model,  # Your XGBoost model
    data = feat,          # Features (predictor variables)
    predict.fun = predict_function
  )
  
  for (var in bin_vars) {
    # Check if the variable has more than one unique value
    if (length(unique(feat[[var]])) > 1) {
      # Calculate ALE for the feature
      ale_effect <- FeatureEffect$new(
        predictor, 
        feature = var,  
        method = "ale"
      )
      
      # Compute net effect from ALE results
      net_effect <- ale_effect$results[2, 2] - ale_effect$results[1, 2]
      ale_val[[var]][i] <- net_effect
    } else {
      # Assign 0 if only one unique value exists
      ale_val[[var]][i] <- 0
    }
  }
}



# Summarize ALE values 
ale_summary <- lapply(ale_val, function(x) {
  data.frame(
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    q2.5 = quantile(x, 0.025, na.rm = TRUE),
    q25 = quantile(x, 0.25, na.rm = TRUE),
    q50 = quantile(x, 0.50, na.rm = TRUE),
    q75 = quantile(x, 0.75, na.rm = TRUE),
    q97.5 = quantile(x, 0.975, na.rm = TRUE)
  )
})

# Convert the list of summaries into a dataframe
ale_summary_df <- do.call(rbind, ale_summary) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Feature") # Create a Feature column from row names

# Print summary dataframe
print(ale_summary_df)

# Create a boxplot
ale_long <- do.call(rbind, lapply(names(ale_val), function(name) {
  data.frame(Feature = name, ALE_Value = ale_val[[name]])
}))

ggplot(ale_long, aes(x = Feature, y = ALE_Value)) +
  geom_boxplot(aes(lower = ..lower.., upper = ..upper.., middle = ..middle.., ymin = ..ymin.., ymax = ..ymax..), 
               stat = "boxplot", width = 0.6, fill = "lightblue") +
  geom_errorbar(aes(ymin = quantile(ALE_Value, 0.025, na.rm = TRUE), 
                    ymax = quantile(ALE_Value, 0.975, na.rm = TRUE)), width = 0.2) +
  theme_minimal() +
  labs(title = "ALE Distribution by Feature", x = "Feature", y = "ALE Effect") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

