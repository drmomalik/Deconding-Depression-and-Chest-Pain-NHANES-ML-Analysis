library(survey)
library(missForest)
library(splitTools)
library(rsample)


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
wd_subset <- wd[, c(responses, predictors, "CDQ001", "DEPR_NEW", "DEPR_NEW", "MEDDEP", "HUQ090", "SDMVPSU", "SDMVSTRA", "MEC15YR", "SEQN")]


# Extract the data from the svrepdesign object
design <- svydesign(
  id = ~SDMVPSU,
  weights = ~MEC15YR,
  strata = ~SDMVSTRA,
  data = wd_subset,
  nest = TRUE
)

survey_design <- subset(design, design$variables["CDQ001"] == 1 & design$variables["DEPR_NEW"] == 1 | design$variables["CDQ001"] == 1 & design$variables["DEPR_NEW"] == 0 & design$variables["MEDDEP"] == 0
                        & design$variables["HUQ090"] == 2 )

strata <- survey_design$variables$SDMVSTRA
clusters <- survey_design$variables$SDMVPSU


survey_design$variables$SDMVSTRA <- NULL
survey_design$variables$SDMVPSU <- NULL
survey_design$variables$SEQN <- NULL
survey_design$variables$CDQ001 <- NULL
survey_design$variables$MEDDEP <- NULL
survey_design$variables$HUQ090 <- NULL


imputed_data <- missForest(xmis=survey_design$variables, verbose = TRUE)
# Extract the completed dataset
wd_subset_imputed <- imputed_data$ximp  # This is the dataset with imputed values
print(imputed_data$OOBerror)

survey_design$variables <- data <- wd_subset_imputed
survey_design$variables$SDMVSTRA <- strata
survey_design$variables$SDMVPSU <- clusters 


#SAVE
both_svy <- survey_design

# create full data replicate samples
rep <- as.svrepdesign(survey_design, type="bootstrap", replicates=100)

repweights <- weights(rep, type="analysis")
repweights <- repweights/mean(repweights) # Stabilize case weights 


########################
library(caret)
library(randomForest)
library(ranger)
library(pROC)
library(MLmetrics)
library(pdp)
library(iml)




data <- survey_design$variables
data <- data[, !duplicated(names(data))]


predictors <- c("AGE_BIN", "RIAGENDR", "RIDRETH1", "SMQ020",
                "SMQ040", "HIQ011", "BPQ020", "BMXBMI", "BPQ080",
                "ALQ130", "PAQMV", "DIDTOT","DMDBORNT",
                "INC3", "DMDEDUC2", "SDDSRVYR", "DUQTOT", "CDQ008",
                "CDQ010", "DEPR_NEW",
                regions)
current_label <- "CAD"
formula <- as.formula(paste(current_label, "~", paste(c(predictors), collapse="+")))

data$CAD <- factor(data$CAD)
levels(data$CAD) <- c("no", "yes")

all_mod <- list()
tune_grid <- expand.grid(
  mtry = c(4, 6, 8),  # Try different values for mtry
  min.node.size = c(1, 5, 10),  # Try different values for min.node.size
  splitrule = "gini"  # Use "gini" for classification (you can also try "extratrees")
)

# Define a custom summary function
customSummary <- function(data, lev = NULL, model = NULL) {
  # 'data$obs' are the observed classes and 'data$pred' are the predicted classes
  sens <- sensitivity(data$pred, data$obs, positive = lev[1])
  spec <- specificity(data$pred, data$obs, negative = lev[2])
  gmean <- sqrt(sens * spec)
  # Return a named vector with your metrics
  out <- c(Sensitivity = sens, Specificity = spec, GMean = gmean)
  return(out)
}

# Use the custom summary function in trainControl:
cv_control <- trainControl(method = "cv", 
                           number = 5, 
                           classProbs = TRUE,
                           summaryFunction = customSummary)


for (rep in 1:ncol(repweights)) {  
  print(paste("On rep sample #:", rep))
  
  # Assign replicate weights
  data$weights <- repweights[, rep]
  
  # Train a Random Forest model with cross-validation and custom tuning grid
  rf_model <- train(formula, 
                    data = data, 
                    method = "ranger", 
                    importance = "permutation",
                    replace = FALSE,
                    weights = data$weights,
                    trControl = cv_control,  
                    num.trees = 100, 
                    tuneGrid = tune_grid,
                    metric = "Specificity"  # Optimize based on your custom 
  )
  all_mod[[rep]] <- rf_model
  
}

# rf_model <- ranger(formula, data=data, importance="permutation", case.weights=data$weights, num.trees = 100, probability = TRUE, class.weights = c(no=1, yes=1000000000000))
# 
# 
# # Initialize vectors (one element per model) to store the metrics
# n_models <- length(all_mod)
# auc_vec   <- numeric(n_models)
# f1_vec    <- numeric(n_models)
# gmean_vec <- numeric(n_models)
# sens_vec  <- numeric(n_models)
# spec_vec  <- numeric(n_models)
# brier_vec <- numeric(n_models)
# 
# for(i in 1:n_models) {
#   cat("On model #:", i, "\n")
#   model <- all_mod[[i]]
#  
#   
#   # Predicted probabilities for the positive class (assumed to be the second column)
#   pred_probs <- predict(model, newdata = filtered_data, type = "prob")[,2]
#   
#   # Predicted classes (for confusion matrix and F1 Score)
#   pred_class <- predict(model, newdata = filtered_data)
#   
#   # The observed classes
#   obs <- filtered_data$CAD
#   
#   # Assuming binary classification: 
#   # Let the positive class be the second level in the factor (e.g., "X1")
#   positive_class <- levels(obs)[2]
#   
#   # Calculate AUC using the predicted probabilities (from pROC package)
#   roc_obj <- roc(obs, pred_probs)
#   auc_vec[i] <- as.numeric(auc(roc_obj))
#   
#   # Compute the confusion matrix using caret
#   cm <- confusionMatrix(as.factor(pred_class), as.factor(obs))
#   
#   sens_vec[i]  <- cm$byClass["Sensitivity"]
#   spec_vec[i]  <- cm$byClass["Specificity"]
#   
#   # Compute F1 Score (using MLmetrics package)
#   f1_vec[i] <- F1_Score(y_pred = pred_class, y_true = obs, positive = positive_class)
#   
#   # Compute Geometric Mean of sensitivity and specificity
#   gmean_vec[i] <- sqrt(as.numeric(cm$byClass["Sensitivity"]) * as.numeric(cm$byClass["Specificity"]))
#   
#   # Compute Brier Score
#   # Convert observed factor to binary (1 if positive class, 0 otherwise)
#   obs_binary <- ifelse(obs == positive_class, 1, 0)
#   brier_vec[i] <- mean((pred_probs - obs_binary)^2)
# }
# 
# # Aggregate the metrics across models (e.g., computing the mean and 95% CI)
# performance <- data.frame(
#   Metric  = c("AUC", "F1", "GMean", "Sensitivity", "Specificity", "Brier"),
#   Mean    = c(mean(auc_vec), mean(f1_vec), mean(gmean_vec), mean(sens_vec), mean(spec_vec), mean(brier_vec)),
#   Lower95 = c(quantile(auc_vec, 0.025), quantile(f1_vec, 0.025),
#               quantile(gmean_vec, 0.025), quantile(sens_vec, 0.025),
#               quantile(spec_vec, 0.025), quantile(brier_vec, 0.025)),
#   Upper95 = c(quantile(auc_vec, 0.975), quantile(f1_vec, 0.975),
#               quantile(gmean_vec, 0.975), quantile(sens_vec, 0.975),
#               quantile(spec_vec, 0.975), quantile(brier_vec, 0.975))
# )
# 
# print(performance)

# Ensure required libraries are loaded
library(pROC)
library(caret)
library(MLmetrics)

# # Split your data into in‐bag and out‐of‐bag sets
# oob_data <- data[data$weights == 0, ]
# ib_data  <- data[data$weights > 0, ]

# Initialize vectors (one element per model) to store the 0.632 bootstrapped metrics
n_models <- length(all_mod)
auc_vec   <- numeric(n_models)
f1_vec    <- numeric(n_models)
gmean_vec <- numeric(n_models)
sens_vec  <- numeric(n_models)
spec_vec  <- numeric(n_models)
brier_vec <- numeric(n_models)

# Define bootstrap weights
w_ib <- 0.368  # weight for in-bag performance
w_oob <- 0.632 # weight for out-of-bag performance

for(i in 1:n_models) {
  cat("Processing model #:", i, "\n")
  model <- all_mod[[i]]
  rep <- repweights[,i]
  
  oob_data <- data[rep == 0, ]
  ib_data <- data[rep > 0, ]
  
  #### Out-of-Bag (OOB) Performance ####
  # Predicted probabilities and classes on OOB data
  pred_probs_oob <- predict(model, newdata = oob_data, type = "prob")[,2]
  pred_class_oob <- predict(model, newdata = oob_data)
  obs_oob <- oob_data$CAD
  
  # Determine positive class (assuming the second level)
  positive_class <- levels(obs_oob)[2]
  
  # AUC for OOB
  roc_obj_oob <- roc(obs_oob, pred_probs_oob)
  auc_oob <- as.numeric(pROC::auc(roc_obj_oob))
  
  # Confusion matrix and related metrics for OOB
  cm_oob <- confusionMatrix(as.factor(pred_class_oob), as.factor(obs_oob))
  sens_oob <- cm_oob$byClass["Sensitivity"]
  spec_oob <- cm_oob$byClass["Specificity"]
  
  # F1 Score for OOB
  f1_oob <- F1_Score(y_pred = pred_class_oob, y_true = obs_oob, positive = positive_class)
  
  # Geometric mean of Sensitivity and Specificity for OOB
  gmean_oob <- sqrt(as.numeric(sens_oob) * as.numeric(spec_oob))
  
  # Brier Score for OOB
  obs_binary_oob <- ifelse(obs_oob == positive_class, 1, 0)
  brier_oob <- mean((pred_probs_oob - obs_binary_oob)^2)
  
  #### In-Bag (IB) Performance ####
  # Predicted probabilities and classes on IB data
  pred_probs_ib <- predict(model, newdata = ib_data, type = "prob")[,2]
  pred_class_ib <- predict(model, newdata = ib_data)
  obs_ib <- ib_data$CAD
  
  # Use the same positive class (assuming levels are consistent)
  positive_class <- levels(obs_ib)[2]
  
  # AUC for IB
  roc_obj_ib <- roc(obs_ib, pred_probs_ib)
  auc_ib <- as.numeric(auc(roc_obj_ib))
  
  # Confusion matrix and related metrics for IB
  cm_ib <- confusionMatrix(as.factor(pred_class_ib), as.factor(obs_ib))
  sens_ib <- cm_ib$byClass["Sensitivity"]
  spec_ib <- cm_ib$byClass["Specificity"]
  
  # F1 Score for IB
  f1_ib <- F1_Score(y_pred = pred_class_ib, y_true = obs_ib, positive = positive_class)
  
  # Geometric mean for IB
  gmean_ib <- sqrt(as.numeric(sens_ib) * as.numeric(spec_ib))
  
  # Brier Score for IB
  obs_binary_ib <- ifelse(obs_ib == positive_class, 1, 0)
  brier_ib <- mean((pred_probs_ib - obs_binary_ib)^2)
  
  #### 0.632 Bootstrap Combination ####
  auc_vec[i]   <- w_ib * auc_ib   + w_oob * auc_oob
  sens_vec[i]  <- w_ib * sens_ib  + w_oob * sens_oob
  spec_vec[i]  <- w_ib * spec_ib  + w_oob * spec_oob
  f1_vec[i]    <- w_ib * f1_ib    + w_oob * f1_oob
  gmean_vec[i] <- w_ib * gmean_ib + w_oob * gmean_oob
  brier_vec[i] <- w_ib * brier_ib + w_oob * brier_oob
}

# Aggregate the bootstrapped metrics across models (for example, computing the mean and 95% CI)
performance <- data.frame(
  Metric  = c("AUC", "F1", "GMean", "Sensitivity", "Specificity", "Brier"),
  Mean    = c(mean(auc_vec), mean(f1_vec), mean(gmean_vec), mean(sens_vec), mean(spec_vec), mean(brier_vec)),
  Lower95 = c(quantile(auc_vec, 0.025), quantile(f1_vec, 0.025),
              quantile(gmean_vec, 0.025), quantile(sens_vec, 0.025),
              quantile(spec_vec, 0.025), quantile(brier_vec, 0.025)),
  Upper95 = c(quantile(auc_vec, 0.975), quantile(f1_vec, 0.975),
              quantile(gmean_vec, 0.975), quantile(sens_vec, 0.975),
              quantile(spec_vec, 0.975), quantile(brier_vec, 0.975))
)

print(performance)


# Create a list to store the importance vectors from each model
importance_list <- lapply(all_mod, function(mod) {
  imp <- mod$finalModel$variable.importance
  # If the importance is a matrix, you might be interested in one measure (e.g., MeanDecreaseGini)
  # Adjust the column name if needed.
  if(is.matrix(imp) && "MeanDecreaseGini" %in% colnames(imp)) {
    imp <- imp[,"MeanDecreaseGini"]
  }
  return(imp)
})

# Combine the importance scores into a matrix:
# Rows = variables, Columns = bootstrap sample (i.e., each model)
importance_mat <- do.call(cbind, importance_list)

# Compute the mean and 95% CI for each variable’s importance score
importance_summary <- data.frame(
  Variable = rownames(importance_mat),
  Mean     = rowMeans(importance_mat, na.rm = TRUE),
  Lower95  = apply(importance_mat, 1, function(x) quantile(x, 0.025, na.rm = TRUE)),
  Upper95  = apply(importance_mat, 1, function(x) quantile(x, 0.975, na.rm = TRUE))
)

importance_summary_sorted <- importance_summary[order(-importance_summary$Mean), ]
print(importance_summary_sorted)

library(ggplot2)

ggplot(importance_summary_sorted, aes(x = reorder(Variable, -Mean), y = Mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = Lower95, ymax = Upper95), width = 0.3) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Variable Importance with 95% CI for predicting CAD in Depression Model",
       x = "Variable", y = "Permuted Performance Difference")



# binary_vars <- c(regions, "CDQ008", "CDQ010")
# 
# # Create a list to store the partial dependence results for each variable.
# # The structure will be: partial_list[[variable]][[replicate]] is the pdp data frame for that replicate.
# partial_list <- list()
# 
# for (var in binary_vars) {
#   partial_list[[var]] <- lapply(seq_along(all_mod), function(i) {
#     # Extract the underlying random forest model from caret's object
#     mod_rf <- all_mod[[i]]
#     # Compute partial dependence for variable 'var'
#     # Setting grid.resolution = 2 ensures that for a binary predictor you get two rows (one per level)
#     pdp_obj <- partial(mod_rf, 
#                        pred.var = var, 
#                        prob = TRUE, 
#                        which.class = 2,   
#                        grid.resolution = 2,
#                        train=data)
#     
#     # Rename the predictor column to "value" for consistency.
#     names(pdp_obj)[names(pdp_obj) == var] <- "value"
#     
#     # Add a column to identify the replicate (optional)
#     pdp_obj$replicate <- i
#     pdp_obj$predictor <- var
#     
#     return(pdp_obj)
#   })
# }
# 
# 
# 
# # Combine all the partial dependence data frames into one
# pdp_combined <- bind_rows(unlist(partial_list, recursive = FALSE))
# 
# # Now, group by the predictor and the level (value) to compute the aggregate statistics.
# agg_partial <- pdp_combined %>%
#   group_by(predictor, value) %>%
#   summarise(
#     Mean    = mean(yhat, na.rm = TRUE),
#     Lower95 = quantile(yhat, 0.025, na.rm = TRUE),
#     Upper95 = quantile(yhat, 0.975, na.rm = TRUE),
#     .groups = "drop"
#   )
# 
# print(agg_partial)
# 
# library(ggplot2)
# 
# ggplot(agg_partial, aes(x = factor(value), y = Mean)) +
#   geom_point(size = 3) +
#   geom_errorbar(aes(ymin = Lower95, ymax = Upper95), width = 0.1) +
#   facet_wrap(~ predictor, scales = "free_y") +
#   labs(x = "Level", y = "Predicted Probability of Outcome",
#        title = "Aggregated Partial Dependence for Binary Predictors") +
#   theme_minimal()
# 
# # For each predictor and replicate, compute the difference:
# # diff = yhat at level "1" minus yhat at level "0"
# diff_df <- pdp_combined %>%
#   group_by(predictor, replicate) %>%
#   summarise(diff = yhat[value == "1"] - yhat[value == "0"]) %>%
#   ungroup()
# diff_summary <- diff_df %>%
#   group_by(predictor) %>%
#   summarise(
#     mean_diff = mean(diff, na.rm = TRUE),
#     lower95  = quantile(diff, 0.025, na.rm = TRUE),
#     upper95  = quantile(diff, 0.975, na.rm = TRUE)
#   ) %>%
#   ungroup()
# 
# print(diff_summary)
# 
# 
# 
# # Sort in descending order
# diff_summary <- diff_summary %>%
#   arrange((mean_diff)) %>%
#   mutate(predictor = factor(predictor, levels = predictor))  # Ensure correct factor order
# 
# # Manually assign names while keeping the original order
# predictor_labels <- c("Right Chest Pain", "Shortness of breath on incline",
#                       "Upper Sternal Pain", "Left Chest Pain", "Epigastric Pain",
#                       "Neck Pain", "Lower Sternal Pain", "Left Arm Pain",
#                       "Severe pain lasting > 30 minutes", "Right Arm Pain")
# 
# # Assign names without disturbing factor order
# diff_summary$predictor <- factor(predictor_labels, levels = rev(predictor_labels))  # Reverse for plotting
# 
# p <- ggplot(diff_summary, aes(x = predictor, y = mean_diff)) +
#   geom_point(size = 3, color = "#2C3E50") +
#   geom_errorbar(aes(ymin = lower95, ymax = upper95), width = 0.2, color = "#2980B9") +
#   geom_text(aes(label = paste("+", sprintf("%.1f%%", mean_diff * 100))), 
#             vjust = -0.5, color = "#2C3E50", size = 5) +
#   coord_flip() +
#   labs(
#     title = "Mean Difference in Predicted CAD Probability",
#     subtitle = "Difference (Level 1 minus Level 0) with 95% Confidence Intervals",
#     x = "Predictor Variable",
#     y = "Mean Difference (Probability)"
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(
#     plot.title = element_text(
#       face = "bold", 
#       margin = ggplot2::margin(t = 0, r = 0, b = 10, l = 0, unit = "pt")
#     ),
#     plot.subtitle = element_text(
#       margin = ggplot2::margin(t = 0, r = 0, b = 15, l = 0, unit = "pt")
#     ),
#     axis.title.y = element_text(
#       margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0, unit = "pt")
#     ),
#     axis.title.x = element_text(
#       margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0, unit = "pt")
#     )
#   )
# 
# # Print the plot
# print(p)
# 


## ALE 
library(iml)
library(dplyr)
library(ggplot2)
library(purrr)   # for unlist(..., recursive = FALSE)

# Define your binary predictor names (make sure they match your data)
binary_vars <- c(regions, "CDQ008", "CDQ010")

# Prediction function
pred_fun <- function(model, newdata) {
  pred <- predict(model, newdata = newdata, type = "prob")[,2]
  
}


# Compute ALE for each variable
ale_list <- list()

for (var in binary_vars) {
  cat("Processing variable:", var, "\n")
  ale_list[[var]] <- lapply(seq_along(all_mod), function(i) {
    # Extract the ranger model from caret
    mod_rf <- all_mod[[i]]
    
    # Create predictor object with corrected model and data
    predictor_obj <- Predictor$new(
      model = mod_rf,
      data = data,  # Original data (no preprocessing applied)
      predict.function = pred_fun
    )
    
    # Compute ALE
    ale_obj <- FeatureEffect$new(
      predictor = predictor_obj,
      feature = var,
      method = "ale",
      grid.size = 2  # For binary variables
    )
    
    # Format results
    ale_data <- ale_obj$results
    names(ale_data)[names(ale_data) == var] <- "value"
    ale_data$replicate <- i
    ale_data$predictor <- var
    
    return(ale_data)
  })
}

# Combine all ALE data frames into one.
ale_combined <- bind_rows(unlist(ale_list, recursive = FALSE))


# Aggregate ALE values for each predictor and level (value).
agg_ale <- ale_combined %>%
  group_by(predictor, value) %>%
  summarise(
    Mean    = mean(.value, na.rm = TRUE),
    Lower95 = quantile(.value, 0.025, na.rm = TRUE),
    Upper95 = quantile(.value, 0.975, na.rm = TRUE),
    .groups = "drop"
  )


# Plot the aggregated ALE values for each binary predictor.
ggplot(agg_ale, aes(x = factor(value), y = Mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Lower95, ymax = Upper95), width = 0.1) +
  facet_wrap(~ predictor, scales = "free_y") +
  labs(x = "Level", y = "Accumulated Local Effect",
       title = "Aggregated ALE for Binary Predictors in Depressed cohort") +
  theme_minimal()
''
# For each predictor and replicate, compute the difference between level "1" and level "0".
# (Assuming the levels are coded as "0" and "1"; adjust if necessary.)
ale_diff_df <- ale_combined %>%
  group_by(predictor, replicate) %>%
  summarise(diff = .value[as.character(value) == "1"] - .value[as.character(value) == "0"]) %>%
  ungroup()

# Summarize the differences to get the mean and 95% confidence intervals.
ale_diff_summary <- ale_diff_df %>%
  group_by(predictor) %>%
  summarise(
    mean_diff = mean(diff, na.rm = TRUE),
    lower95  = quantile(diff, 0.025, na.rm = TRUE),
    upper95  = quantile(diff, 0.975, na.rm = TRUE)
  ) %>%
  ungroup()

# First, sort and set the factor order as before:
ale_diff_summary <- ale_diff_summary %>%
  arrange(mean_diff) %>%
  mutate(predictor = factor(predictor, levels = predictor))

# Create a named vector with the new labels:
new_labels <- c(
  "CDQ009B" = "Right Chest Pain",
  "CDQ010"  = "Shortness of breath on incline",
  "CDQ009D" = "Upper Sternal Pain",
  "CDQ009F" = "Left Chest Pain",
  "CDQ009H" = "Epigastric Pain",
  "CDQ009C" = "Neck Pain",
  "CDQ009E" = "Lower Sternal Pain",
  "CDQ009G" = "Left Arm Pain",
  "CDQ008"  = "Severe pain lasting > 30 minutes",
  "CDQ009A" = "Right Arm Pain"
)

# Now, update the levels of the factor directly:
current_levels <- levels(ale_diff_summary$predictor)
levels(ale_diff_summary$predictor) <- new_labels[current_levels]


# Define the mapping for outcome labels
new_labels <- c(
  "CDQ009D" = "Upper Sternal Pain",
  "CDQ009F" = "Left Chest Pain",
  "CDQ009H" = "Epigastric Pain",
  "CDQ009C" = "Neck Pain",
  "CDQ009E" = "Lower Sternal Pain",
  "CDQ009G" = "Left Arm Pain",
  "CDQ009B"  = "Right Chest Pain",
  "CDQ009A" = "Right Arm Pain"
)

# Create the plot
p <- ggplot(plot_df, aes(x = Outcome, fill = DEPR_NEW)) +
  # Whiskers: thin error bars representing ±2 SD from the mean
  geom_errorbar(aes(ymin = LowerWhisker, ymax = UpperWhisker), 
                width = 0.2, 
                position = position_dodge(0.5)) +
  # Box: thick error bar representing ±1 SD from the mean
  geom_errorbar(aes(ymin = LowerBox, ymax = UpperBox), 
                width = 0.5, 
                position = position_dodge(0.5), 
                size = 2) +
  # Mean as a point marker
  geom_point(aes(y = Mean), 
             position = position_dodge(0.5), 
             size = 3, shape = 21, color = "black") +
  labs(x = "Outcome",
       y = "Proportion",
       fill = "DEPR_NEW",
       title = "Chest Pain Location Survey-Weighted Proportions by Depression") +
  # Use the new labels for x-axis values
  scale_x_discrete(labels = new_labels) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, family = "Times New Roman"),
    text = element_text(family = "Times New Roman")
  )

print(p)


##### Stratified ALE by DEPR_NEW

## ALE with stratification by DEPR_NEW
library(iml)
library(dplyr)
library(ggplot2)
library(purrr)   # for unlist(..., recursive = FALSE)

# Define your binary predictor names (make sure they match your data)
binary_vars <- c(regions, "CDQ008", "CDQ010")

# Prediction function
pred_fun <- function(model, newdata) {
  predict(model, newdata = newdata, type = "prob")[,2]
}

# Compute ALE for each variable stratified by depression status ("DEPR_NEW")
ale_list <- list()

for (var in binary_vars) {
  cat("Processing variable:", var, "\n")
  # For each predictor, run ALE over all model replicates
  ale_list[[var]] <- lapply(seq_along(all_mod), function(i) {
    mod_rf <- all_mod[[i]]
    
    # Loop over depression groups (e.g., DEPR_NEW == 0 and DEPR_NEW == 1)
    depression_levels <- unique(data$DEPR_NEW)
    group_results <- lapply(depression_levels, function(dep_val) {
      # Subset the data to the current depression group
      data_subset <- data %>% filter(DEPR_NEW == dep_val)
      
      # Create iml predictor on the subset
      predictor_obj <- Predictor$new(
        model = mod_rf,
        data = data,
        predict.function = pred_fun
      )
      
      # Compute ALE for the current predictor variable (grid.size = 2 works for binary variables)
      ale_obj <- FeatureEffect$new(
        predictor = predictor_obj,
        feature = c(var, "DEPR_NEW"),
        grid.size = 2
      )
      
      ale<-ALEPlot(data[,c(24,29)], mod_rf, pred.fun=pred_fun, J=, K=20, NA.plot=TRUE)
      # Format the ALE output
      ale_data <- ale_obj$results
      names(ale_data)[names(ale_data) == var] <- "value"
      ale_data$replicate <- i
      ale_data$predictor <- var
      ale_data$depression <- as.character(dep_val)  # keep as character or convert to factor as needed
      
      return(ale_data)
    })
    # Combine the ALE results from both depression groups for this replicate
    bind_rows(group_results)
  })
}

# Combine all ALE data frames into one.
ale_combined <- bind_rows(unlist(ale_list, recursive = FALSE))

# Aggregate ALE values for each predictor, depression group, and level (value).
agg_ale <- ale_combined %>%
  group_by(predictor, depression, value) %>%
  summarise(
    Mean    = mean(.value, na.rm = TRUE),
    Lower95 = quantile(.value, 0.025, na.rm = TRUE),
    Upper95 = quantile(.value, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# Plot the aggregated ALE values for each binary predictor stratified by depression.
ggplot(agg_ale, aes(x = factor(value), y = Mean, color = depression, group = depression)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Lower95, ymax = Upper95), width = 0.1) +
  facet_wrap(~ predictor, scales = "free_y") +
  labs(x = "Level", y = "Accumulated Local Effect",
       title = "Aggregated ALE for Binary Predictors Stratified by Depression (DEPR_NEW)") +
  theme_minimal()

# --- Compute ALE Differences ---

# For each predictor, replicate, and depression group, compute the difference
# between ALE at level "1" and level "0" (adjust levels if necessary)
ale_diff_df <- ale_combined %>%
  group_by(predictor, replicate, depression) %>%
  summarise(diff = .value[as.character(value) == "1"] - .value[as.character(value) == "0"]) %>%
  ungroup()

# Summarize the differences to get the mean and 95% confidence intervals for each predictor and depression group.
ale_diff_summary <- ale_diff_df %>%
  group_by(predictor, depression) %>%
  summarise(
    mean_diff = mean(diff, na.rm = TRUE),
    lower95  = quantile(diff, 0.025, na.rm = TRUE),
    upper95  = quantile(diff, 0.975, na.rm = TRUE)
  ) %>%
  ungroup()

ale_diff_summary <- ale_diff_summary %>%
  arrange(mean_diff) %>%
  mutate(predictor = factor(predictor, levels = unique(predictor)))


# Create a named vector with the new labels (update these labels as needed)
new_labels <- c(
  "CDQ009B" = "Right Chest Pain",
  "CDQ010"  = "Shortness of breath on incline",
  "CDQ009D" = "Upper Sternal Pain",
  "CDQ009F" = "Left Chest Pain",
  "CDQ009H" = "Epigastric Pain",
  "CDQ009C" = "Neck Pain",
  "CDQ009E" = "Lower Sternal Pain",
  "CDQ009G" = "Left Arm Pain",
  "CDQ008"  = "Severe pain lasting > 30 minutes",
  "CDQ009A" = "Right Arm Pain"
)

current_levels <- levels(ale_diff_summary$predictor)
levels(ale_diff_summary$predictor) <- new_labels[current_levels]

# Create the plot of the mean differences (ALE: Level 1 minus Level 0) with 95% CIs,
# faceted by depression status.
p <- ggplot(ale_diff_summary, aes(x = predictor, y = mean_diff)) +
  geom_point(size = 3, color = "#2C3E50") +
  geom_errorbar(aes(ymin = lower95, ymax = upper95), width = 0.2, color = "#2980B9") +
  geom_text(aes(label = paste("+", sprintf("%.1f%%", mean_diff * 100))),
            vjust = -0.5, color = "#2C3E50", size = 5) +
  facet_wrap(~ depression, scales = "free_y") +
  coord_flip() +
  labs(
    title = "Mean Difference in Predicted CAD Probability (ALE) Stratified by Depression",
    subtitle = "Difference (Level 1 minus Level 0) with 95% Confidence Intervals",
    x = "Predictor Variable",
    y = "Mean Difference (Probability)"
  ) +
  theme_minimal(base_size = 14) 
 

# Print the plot.
print(p)


#####


library(flextable)

# Create the flextable
performance_flex <- flextable(performance) %>%
  set_header_labels(Metric = "Metric", Mean = "Mean (95% CI)", Lower95 = "95% CI Lower", Upper95 = "95% CI Upper") %>%
  compose(j = "Metric", value = as_paragraph((Metric))) %>%
  compose(j = "Mean", value = as_paragraph(format(Mean, digits = 3))) %>%
  compose(j = "Lower95", value = as_paragraph(format(Lower95, digits = 3))) %>%
  compose(j = "Upper95", value = as_paragraph(format(Upper95, digits = 3))) %>%
  bold(i = 1, part = "header") %>%
  align(align = "center", part = "all") %>%
  autofit()

# Display the flextable
performance_flex
