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

split <- initial_split(survey_design$variables, strata = "CAD", prop = 0.75)
train_indices <- split$in_id
test_indices <- setdiff(1:nrow(split$data), train_indices)
train_design <- survey_design[train_indices,]
test_design <- survey_design[test_indices,]

train_data <- train_design$variables
test_data <- test_design$variables

# Create a train rep and test rep designs for validation of ECC model 
train_repwt <- as.svrepdesign(train_design, type="bootstrap", replicates=10)
test_rep <- as.svrepdesign(test_design, type="bootstrap", replicates=10)

# Need to create actual weights in order to bootstrap while training (do not need to do for testing)
train_repweights <- weights(train_repwt, type = "analysis") / mean(weights(train_repwt, type = "analysis"))


###################################################

current_label <- "CAD"
formula <- as.formula(paste(current_label, "~", paste(c(predictors), collapse="+")))
data <- train_data
data[[current_label]] <- as.factor(data[[current_label]])

data$weights <- train_repweights[,2]

data_minority <- data %>% filter(CAD == 1)
data_majority <- data %>% filter(CAD == 0)
ratio <- nrow(data_majority)/nrow(data_minority)

beta <- 2
majority_sampled <- data_majority %>% sample_n(nrow(data_minority), replace = TRUE, weight=weights)
minority_sampled <- data_minority %>% sample_n(nrow(data_minority), replace = TRUE, weight=weights)
balanced_data <- bind_rows(minority_sampled, majority_sampled)
balanced_data[[current_label]] <- as.factor(balanced_data[[current_label]])

model <- imbalanced(formula, data=data, case.wt = data$weights, 
                    method="rfq", splitrule = "auc", ntree =3000, 
                    nodesize=5, fast=TRUE, forest=TRUE, importance="permute",
                    samptype="swr", perf.type = "gmean")

pred<-predict(model, newdata=test_data, outcome="test")
print(model)
print(pred)
th <- get.imbalanced.optimize(model)["threshold"]
get.imbalanced.performance(model, threshold= th)
get.imbalanced.performance(pred, threshold = th)

synth <- synthetic.rfsrc(formula, data=data, case.wt=data$weights, ntree=2000, 
                         newdata = test_data, splitrule="auc")


mod <- svyglm(formula, design=train_design, family=quasibinomial())
predict <- predict(mod, newdata=test_data, type="response")
predictions <- ifelse(as.numeric(predict) >= 0.5, 1, 0)
true_label <- test_data$CAD

auc <- auc(true_label, predictions)

find.interaction(model, verbose=TRUE)
vimp <- vimp(model, importance="permute", joint=FALSE)
vimp_ci <- subsample.rfsrc(model, B=50, verbose=TRUE, importance="permute", joint=FALSE)
par(mfrow=c(1,1))
plot.subsample(vimp_ci, alpha=0.05)
print(model)


plot.variable(model, c(predictors), partial=TRUE, sorted=TRUE)



####################

current_label <- "CAD"
formula <- as.formula(paste(current_label, "~", paste(c(predictors), collapse="+")))
# Initialize an empty randomForest object to store all trees
all <- NULL

# Loop to train 20 models and combine them
for (i in 1:ncol(train_repweights)) {
  print(paste("On replicate #:", i))
  train_data$weights <- train_repweights[,i]+0.001
  # Train a randomForest model with case weights
  model <- randomForest(
    formula, 
    data = train_data, 
    weights = train_data$weights,
    ntree = 100,
    nodesize=5
  )
  
  # Combine models
  if (is.null(all)) {
    all <- model  # Initialize with the first model
  } else {
    all <- combine(all, model)  # Combine subsequent models
  }
}

# Print the final combined model
print(all)

true_label <- test_data$CAD
predictions <- predict(all, newdata=test_data, type="response")

predictions <- as.numeric(as.character(predictions))
library(pROC)

auc <- auc(true_label, predictions)

print(auc)

# Calculate F1 Score
calculate_f1 <- function(true_labels, predicted_labels) {
  # Confusion matrix
  tp <- sum(true_labels == 1 & predicted_labels == 1)  # True positives
  fp <- sum(true_labels == 0 & predicted_labels == 1)  # False positives
  fn <- sum(true_labels == 1 & predicted_labels == 0)  # False negatives
  
  # Precision and Recall
  precision <- tp / (tp + fp)
  recall <- tp / (tp + fn)
  
  # F1 Score
  f1 <- 2 * (precision * recall) / (precision + recall)
  return(data.frame(f1=f1, tp=tp, fp=fp, fn=fn, precision=precision, recall=recall))
}

f1_score <- calculate_f1(true_label, predictions)
print(f1_score)






########################
library(caret)
library(randomForest)
library(pROC)
library(MLmetrics)
library(pdp)



data <- survey_design$variables


predictors <- c("AGE_BIN", "RIAGENDR", "RIDRETH1", "SMQ020",
                "SMQ040", "HIQ011", "BPQ020", "BMXBMI", "BPQ080",
                "ALQ130", "PAQMV", "DIDTOT","DMDBORNT",
                "INC3", "DMDEDUC2", "SDDSRVYR", "DUQTOT", "CDQ008",
                "CDQ010",
                regions)
current_label <- "CAD"
formula <- as.formula(paste(current_label, "~", paste(c(predictors), collapse="+")))

all_model <- list()
for (rep in 1:ncol(repweights)) {  
  print(paste("On rep sample #:", rep))
  # Define cross-validation method
  cv_control <- trainControl(method = "cv", number = 5)  # 5-fold cross-validation
  data$weights <- repweights[,rep]
  # Train a Random Forest model with cross-validation
  rf_model <- train(formula, 
                    data = data, 
                    method = "rf", 
                    importance = TRUE,
                    weights = data$weights,
                    trControl = cv_control,
                    ntree=50,
                    nodesize=10,
                    tuneLength = 3)  # Tuning grid size
  
  all_model[[rep]] <- rf_model
}


# Initialize vectors (one element per model) to store the metrics
n_models <- length(all_model)
auc_vec   <- numeric(n_models)
f1_vec    <- numeric(n_models)
gmean_vec <- numeric(n_models)
sens_vec  <- numeric(n_models)
spec_vec  <- numeric(n_models)
brier_vec <- numeric(n_models)

for(i in 1:n_models) {
  print(paste("On model #:", i))
  model <- all_model[[i]]
  
  # Predicted class
  preds <- predict(model)
  preds <- as.numeric(as.character(preds))
  
  # The observed classes
  obs <- data$CAD
  
  # Assuming binary classification: 
  # Let the positive class be the second level in the factor
  positive_class <- levels(obs)[2]
  
  
  # Calculate AUC using the pROC package
  roc_obj <- roc(obs, preds)
  auc_vec[i] <- as.numeric(auc(roc_obj))
  
  # Compute the confusion matrix using caret
  cm <- confusionMatrix(as.factor(preds), as.factor(obs))
  
  sens_vec[i]  <- cm$byClass["Sensitivity"]
  spec_vec[i]  <- cm$byClass["Specificity"]
  
  # Compute F1 Score (using MLmetrics)
  f1_vec[i] <- F1_Score(y_pred = preds, y_true = obs, positive = positive_class)
  
  # Compute Geometric Mean of sensitivity and specificity
  gmean_vec[i] <- sqrt(as.numeric(cm$byClass["Sensitivity"]) * as.numeric(cm$byClass["Specificity"]))
  
  # Compute Brier Score
  # Convert observed factor to binary (0/1)
  obs_binary <- ifelse(obs == positive_class, 1, 0)
  brier_vec[i] <- mean((preds - obs_binary)^2)
}

# Now aggregate the metrics across models. For example, computing the mean and 95% CI:
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
importance_list <- lapply(all_model, function(mod) {
  imp <- mod$finalModel$importance
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
  labs(title = "Variable Importance with 95% CI",
       x = "Variable", y = "Mean Decrease in Gini")



binary_vars <- c(regions, "CDQ008", "CDQ010")

# Create a list to store the partial dependence results for each variable.
# The structure will be: partial_list[[variable]][[replicate]] is the pdp data frame for that replicate.
partial_list <- list()

for (var in binary_vars) {
  partial_list[[var]] <- lapply(seq_along(all_model), function(i) {
    # Extract the underlying random forest model from caret's object
    mod_rf <- all_model[[i]]
    # Compute partial dependence for variable 'var'
    # Setting grid.resolution = 2 ensures that for a binary predictor you get two rows (one per level)
    pdp_obj <- partial(mod_rf, 
                       pred.var = var, 
                       prob = TRUE, 
                       which.class = 2,   
                       grid.resolution = 2,
                       train=data)
    
    # Rename the predictor column to "value" for consistency.
    names(pdp_obj)[names(pdp_obj) == var] <- "value"
    
    # Add a column to identify the replicate (optional)
    pdp_obj$replicate <- i
    pdp_obj$predictor <- var
    
    return(pdp_obj)
  })
}



# Combine all the partial dependence data frames into one
pdp_combined <- bind_rows(unlist(partial_list, recursive = FALSE))

# Now, group by the predictor and the level (value) to compute the aggregate statistics.
agg_partial <- pdp_combined %>%
  group_by(predictor, value) %>%
  summarise(
    Mean    = mean(yhat, na.rm = TRUE),
    Lower95 = quantile(yhat, 0.025, na.rm = TRUE),
    Upper95 = quantile(yhat, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

print(agg_partial)

library(ggplot2)

ggplot(agg_partial, aes(x = factor(value), y = Mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Lower95, ymax = Upper95), width = 0.1) +
  facet_wrap(~ predictor, scales = "free_y") +
  labs(x = "Level", y = "Predicted Probability of Outcome",
       title = "Aggregated Partial Dependence for Binary Predictors") +
  theme_minimal()

# For each predictor and replicate, compute the difference:
# diff = yhat at level "1" minus yhat at level "0"
diff_df <- pdp_combined %>%
  group_by(predictor, replicate) %>%
  summarise(diff = yhat[value == "1"] - yhat[value == "0"]) %>%
  ungroup()
diff_summary <- diff_df %>%
  group_by(predictor) %>%
  summarise(
    mean_diff = mean(diff, na.rm = TRUE),
    lower95  = quantile(diff, 0.025, na.rm = TRUE),
    upper95  = quantile(diff, 0.975, na.rm = TRUE)
  ) %>%
  ungroup()

print(diff_summary)



# Sort in descending order
diff_summary <- diff_summary %>%
  arrange((mean_diff)) %>%
  mutate(predictor = factor(predictor, levels = predictor))  # Ensure correct factor order

# Manually assign names while keeping the original order
predictor_labels <- c("Epigastric Pain", "Right Arm Pain",
                      "Shortness of breath on incline", "Neck Pain", "Lower Sternal Pain",
                      "Left Arm Pain", "Right Chest Pain", "Left Chest Pain",
                      "Upper Sternal Pain", "Severe Pain in lasting > 30 minutes")

# Assign names without disturbing factor order
diff_summary$predictor <- factor(predictor_labels, levels = rev(predictor_labels))  # Reverse for plotting

p <- ggplot(diff_summary, aes(x = predictor, y = mean_diff)) +
  geom_point(size = 3, color = "#2C3E50") +
  geom_errorbar(aes(ymin = lower95, ymax = upper95), width = 0.2, color = "#2980B9") +
  geom_text(aes(label = paste("+", sprintf("%.1f%%", mean_diff * 100))), 
            vjust = -0.5, color = "#2C3E50", size = 5) +
  coord_flip() +
  labs(
    title = "Mean Difference in Predicted CAD Probability without History of Depression",
    subtitle = "Difference (Level 1 minus Level 0) with 95% Confidence Intervals",
    x = "Predictor Variable",
    y = "Mean Difference (Probability)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(
      face = "bold", 
      margin = ggplot2::margin(t = 0, r = 0, b = 10, l = 0, unit = "pt")
    ),
    plot.subtitle = element_text(
      margin = ggplot2::margin(t = 0, r = 0, b = 15, l = 0, unit = "pt")
    ),
    axis.title.y = element_text(
      margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0, unit = "pt")
    ),
    axis.title.x = element_text(
      margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0, unit = "pt")
    )
  )

# Print the plot
print(p)






# 
# # View results
# print(rf_model1)
# plot(rf_model1)  # Visualize performance across different hyperparameters
# 
# pred<-predict(rf_model1, newdata=data)
# predictions<-as.numeric(as.character(pred))
# predictions <- as.numeric(as.character(combo$predicted))
# true_label <- data$CAD
# 
# library(pROC)
# 
# auc <- auc(true_label, predictions)
# 
# print(auc)
# 
# # Calculate F1 Score
# calculate_f1 <- function(true_labels, predicted_labels) {
#   # Confusion matrix
#   tp <- sum(true_labels == 1 & predicted_labels == 1)  # True positives
#   fp <- sum(true_labels == 0 & predicted_labels == 1)  # False positives
#   fn <- sum(true_labels == 1 & predicted_labels == 0)  # False negatives
#   
#   # Precision and Recall
#   precision <- tp / (tp + fp)
#   recall <- tp / (tp + fn)
#   
#   # F1 Score
#   f1 <- 2 * (precision * recall) / (precision + recall)
#   return(data.frame(f1=f1, tp=tp, fp=fp, fn=fn, precision=precision, recall=recall))
# }
# 
# f1_score <- calculate_f1(true_label, predictions)
# print(f1_score)
# 
# 
