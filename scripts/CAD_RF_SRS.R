library(survey)
library(tableone)  # For balance checking
library(MatchIt)   # For propensity score estimation
library(splines)
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

predictors <- c("AGE_BIN","AGE_BIN*RIAGENDR", "RIAGENDR", "RIDRETH1", "SMQ020",
                "SMQ040", "HIQ011", "BPQ020", "log(BMXBMI)", "BPQ080",
                "log(ALQ130)", "PAQMV", "DIDTOT","DMDBORNT",
                "INC3", "DMDEDUC2", "SDDSRVYR", "DUQTOT")
covariates <- c("AGE_BIN", "RIAGENDR", "RIDRETH1", "SMQ020",
                "SMQ040", "HIQ011", "BPQ020", "BMXBMI", "BPQ080",
                "ALQ130", "PAQMV", "DIDTOT","DMDBORNT",
                "INC3", "DMDEDUC2", "SDDSRVYR", "DUQTOT","ps")
current_label <- "DEPR_NEW"
formula <- as.formula(paste(current_label, "~", paste(c(predictors), collapse="+")))

data <- survey_design$variables
# Step 2: Estimate propensity scores with survey weights
ps_model <- svyglm(formula, design = survey_design, family = quasibinomial())
data$ps <- predict(ps_model, type = "response")

# Step 3: Match within survey strata (do not need to use weights)
matched_data <- matchit(
  DEPR_NEW ~ ps,
  data = data,
  method = "nearest", 
  exact = ~SDMVSTRA  # Match within strata
)

# Remove rows where SDMVSTRA is in the specified strata
matched_df <- match.data(matched_data) %>%
  filter(!SDMVSTRA %in% c(51, 55, 63, 57, 70, 74, 87))
levels(matched_df$CAD) <- c("no", "yes")


# Assign treated units' weights to their matched controls (transfer weights)
matched_df$weight_transferred <- ifelse(
  matched_df$DEPR_NEW == 1, 
  matched_df$MEC15YR, 
  matched_df$MEC15YR[matched_df$subclass]  # Transfer within subclass
)


# Step 4: Analyze matched data with survey weights
matched_design <- svydesign(
  ids = ~SDMVPSU, 
  strata = ~SDMVSTRA, 
  weights = ~MEC15YR, 
  data = matched_df,
  nest=TRUE
)




original_design <- survey_design
original_data <- data

library(survey)
library(ggplot2)
library(cobalt)

# Check matched covariate balance 
check_survey_balance <- function(design, treatment_var, covariates) {
  
  # Convert design to dataframe with weights
  df <- model.frame(design)
  
  # Ensure numeric covariates only
  numeric_covars <- covariates[sapply(df[covariates], is.numeric)]
  
  balance_table <- lapply(numeric_covars, function(x) {
    # Get survey means
    means <- svyby(as.formula(paste("~", x)),
                   as.formula(paste("~", treatment_var)),
                   design, svymean)
    
    # Get survey variances
    vars <- svyby(as.formula(paste("~", x)),
                  as.formula(paste("~", treatment_var)),
                  design, svyvar)
    
    # Calculate SMD
    smd <- (means[2, x] - means[1, x]) / sqrt((vars[1, x] + vars[2, x])/2)
    
    data.frame(
      Variable = x,
      SMD = smd,
      Treatment_Mean = means[2, x],
      Control_Mean = means[1, x]
    )
  }) |> dplyr::bind_rows()
  
  # Handle factor variables separately
  factor_covars <- covariates[sapply(df[covariates], is.factor)]
  if(length(factor_covars) > 0) {
    factor_balance <- lapply(factor_covars, function(f) {
      lvls <- levels(df[[f]])
      lapply(lvls, function(l) {
        df$temp <- as.numeric(df[[f]] == l)
        design_temp <- update(design, temp = df$temp)
        means <- svyby(~temp, as.formula(paste("~", treatment_var)),
                       design_temp, svymean)
        data.frame(
          Variable = paste0(f, ":", l),
          SMD = (means[2, "temp"] - means[1, "temp"]) /
            sqrt((means[1, "temp"]*(1-means[1, "temp"]) +
                    means[2, "temp"]*(1-means[2, "temp"]))/2),
          Treatment_Mean = means[2, "temp"],
          Control_Mean = means[1, "temp"]
        )
      }) |> dplyr::bind_rows()
    }) |> dplyr::bind_rows()
    
    balance_table <- dplyr::bind_rows(balance_table, factor_balance)
  }
  
  return(balance_table)
}

# 3. Create balance tables

original_balance <- check_survey_balance(original_design, "DEPR_NEW", covariates)
matched_balance <- check_survey_balance(matched_design, "DEPR_NEW", covariates)

# 4. Combine balance metrics
final_balance <- merge(original_balance, matched_balance, by = "Variable",
                       suffixes = c("_before", "_after"))
final_balance <- filter(final_balance, !Variable == "ps")

# 5. Visualize balance using Love plot
library(ggplot2)
library(ggrepel)  # Load ggrepel for better label placement

love.plot <- ggplot(final_balance, aes(x = SMD_before, y = SMD_after)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_hline(yintercept = c(-0.1, 0.1), color = "red") +
  geom_text_repel(aes(label = Variable), size = 3, family = "Times New Roman", max.overlaps = 15) +  # Use ggrepel to avoid overlap
  ggtitle("Love Plot of Standardized Mean Differences") +
  xlab("SMD Before Matching") + 
  ylab("SMD After Matching") +
  theme_minimal(base_family = "Times New Roman") +  
  theme(
    text = element_text(family = "Times New Roman", size = 16),  
    axis.title.x = element_text(size = 18, face = "bold"),  
    axis.title.y = element_text(size = 18, face = "bold"),  
    axis.text = element_text(size = 16),  
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)  
  )

love.plot



# # 6. Distribution of propensity scores
# ps_plot <- ggplot(matched_data, aes(x = ps, fill = factor(treatment))) +
#   geom_density(alpha = 0.5) +
#   ggtitle("Propensity Score Distribution After Matching") +
#   xlab("Propensity Score") + ylab("Density") +
#   theme_minimal()
# 
# # 7. Balance metrics table
# balance_metrics <- final_balance
# balance_metrics$Improvement <- balance_metrics$SMD_before - balance_metrics$SMD_after
# 
# # 8. Variance ratios (important for balance assessment)
# var_ratio <- function(var) {
#   svyvar(~var, subset(original_design, treatment == 1)) /
#     svyvar(~var, subset(original_design, treatment == 0))
# }
# 
# balance_metrics$VarRatio_before <- sapply(covariates, var_ratio)
# balance_metrics$VarRatio_after <- sapply(covariates, 
#                                          function(x) svyvar(~x, matched_design) )
# 
# # Print results
# print(list(
#   balance_table = balance_metrics,
#   love_plot = love.plot,
#   ps_plot = ps_plot
# ))

### Make replicate weights and datasets for analysis

rep <- as.svrepdesign(matched_design, type="bootstrap", replicates=250)
repweights <- weights(rep, type="analysis")
repweights <- repweights/mean(repweights) # Stabilize case weights 



# Create dataset with each repweight
all_data <- list()
for (i in 1:ncol(repweights)) {
  df <- matched_df
  df$weights <- repweights[,i]
  all_data[[i]] <- df
}

### Run meta-model for depression cohort on CAD
vars <- c("AGE_BIN", "RIAGENDR", "RIDRETH1", "SMQ020",
          "SMQ040", "HIQ011", "BPQ020", "BMXBMI", "BPQ080",
          "ALQ130", "PAQMV", "DIDTOT","DMDBORNT",
          "INC3", "DMDEDUC2", "SDDSRVYR", "DUQTOT", "CDQ008",
          "CDQ010",
          regions)
outcome <- "CAD"
formula_new <- as.formula(paste(outcome, "~", paste(c(vars), collapse="+")))


all_m_nw <- list()
tune_grid <- expand.grid(
  mtry = c(8, 16),  # Try different values for mtry
  min.node.size = c(1,5,10),  # Try different values for min.node.size4
  splitrule = "gini"  # Use "gini" for classification (you can also try "extratrees")
)

# Define a custom summary function
customSummary <- function(data, lev = NULL, model = NULL) {
  # 'data$obs' are the observed classes and 'data$pred' are the predicted classes
  sens <- caret::sensitivity(data$pred, data$obs, positive = lev[1])
  spec <- caret::specificity(data$pred, data$obs, negative = lev[2])
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
  tempdata <- all_data[[rep]]
  newdata <- tempdata[tempdata$DEPR_NEW==1,]
  
  # Train a Random Forest model with cross-validation and custom tuning grid
  rf_model <- caret::train(formula_new, 
                           data = newdata, 
                           method = "ranger", 
                           importance = "permutation",
                           replace = TRUE,
                           trControl = cv_control, 
                           num.trees = 100, 
                           tuneGrid = tune_grid,
                           metric = "Specificity"  # Optimize based on your custom 
  )
  all_m_nw[[rep]] <- rf_model
  
}

# Initialize vectors (one element per model) to store the 0.632 bootstrapped metrics
n_models <- length(all_m_nw)
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
  model <- all_m_nw[[i]]
  rep <- repweights[,i]
  
  tempdata <- all_data[[i]]
  newdata <- tempdata[tempdata$DEPR_NEW==1,]
  new_dat<-newdata
  levels(new_dat)<-c("0","1")
  oob_data <- newdata[newdata$weights == 0, ]
  ib_data <- newdata[newdata$weights > 0, ]
  
  #### Out-of-Bag (OOB) Performance ####
  # Predicted probabilities and classes on OOB data
  pred_probs_oob <- predict(model, newdata = oob_data, type = "prob")[,2]
  pred_class_oob <- predict(model, newdata = oob_data)
  levels(pred_class_oob)<-c(0,1)
  obs_oob <- oob_data$CAD
  levels(obs_oob)<-c(0,1)
  
  # Determine positive class (assuming the second level)
  positive_class <- levels(obs_oob)[2]
  
  # AUC for OOB
  roc_obj_oob <- pROC::roc(obs_oob, pred_probs_oob)
  auc_oob <- as.numeric(pROC::auc(roc_obj_oob))
  
  # Confusion matrix and related metrics for OOB
  cm_oob <- table(pred_class_oob, obs_oob)
  sens_oob <- cm_oob["1", "1"] / sum(cm_oob["1", ])  # TP / (TP + FN)
  spec_oob <- cm_oob["0", "0"] / sum(cm_oob["0", ])  # TN / (TN + FP)
  
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
  levels(pred_class_ib)<-c(0,1)
  obs_ib <- ib_data$CAD
  levels(obs_ib)<-c(0,1)
  
  # Use the same positive class (assuming levels are consistent)
  positive_class <- levels(obs_ib)[2]
  
  # AUC for IB
  roc_obj_ib <- pROC::roc(obs_ib, pred_probs_ib)
  auc_ib <- as.numeric(pROC::auc(roc_obj_ib))
  
  # Confusion matrix and related metrics for IB
  cm_ib <- table(pred_class_ib, obs_ib)
  sens_ib <- cm_ib["1", "1"] / sum(cm_ib["1", ])  # TP / (TP + FN)
  spec_ib <- cm_ib["0", "0"] / sum(cm_ib["0", ])  # TN / (TN + FP)
  
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


### VIMP

# Create a list to store the importance vectors from each model
importance_list <- lapply(all_m_nw, function(mod) {
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


#### Sig VIMP
# Filter to keep only variables where the 95% CI does not cross zero
significant_vimp <- importance_summary_sorted[importance_summary_sorted$Lower95 > 0, ]

labels0 <- c("CDQ009G1"="Left Arm Pain",
             "INC32"="Income as proportion of poverty threshold: 1-1.99",
             "BMXBMI"="BMI",
             "AGE_BIN3"="Age: 60-69 years",
             "INC33"="Income as proportion of poverty threshold: > 2",
             "DUQTOT1"="Use of illicit drugs",
             "ALQ130"="Alcohol use (drinks/day)",
             "SMQ0403"="Current smoker",
             "RIAGENDR2"="Male",
             "DIDTOT4"="Diabetes: Insulin",
             "AGE_BIN5"="Age: > 80 years",
             "BPQ0202"="History of Hypertension",
             "BPQ0802"="History of Dyslipidemia",
             "CDQ0081"="Severe chest pain for > 30 minutes")


ggplot(significant_vimp, aes(x = reorder(Variable, -Mean), y = Mean)) +
  geom_point(size = 3, color = "red2") +  # Make points larger and color them blue
  geom_errorbar(aes(ymin = Lower95, ymax = Upper95), 
                width = 0.3, 
                size = 1, 
                color = "red2") +  # Thicker error bars in dark red
  coord_flip() +
  theme_minimal(base_family = "Times New Roman", base_size = 16) +  # Set Times New Roman globally
  labs(
    title = "Significant Variable Importance with 95% CI \n for Predicting CAD in Depression Model",
    x = "Variable", 
    y = "Permuted Performance Difference"
  ) +
  scale_x_discrete(labels = labels0) +  # Ensure labels are properly displayed
  theme(
    axis.text.x = element_text(size = 16, family = "Times New Roman"),  
    axis.text.y = element_text(size = 16, family = "Times New Roman"),  
    axis.title.x = element_text(size = 18, face = "bold", family = "Times New Roman"),  
    axis.title.y = element_text(size = 18, face = "bold", family = "Times New Roman"),  
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5, family = "Times New Roman")  
  )



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
  ale_list[[var]] <- lapply(seq_along(all_m_nw), function(i) {
    # Extract the ranger model from caret
    mod_rf <- all_m_nw[[i]]
    tempdata <- all_data[[i]]
    newdata <- tempdata[tempdata$DEPR_NEW==1,]
    
    # Create predictor object with corrected model and data
    predictor_obj <- Predictor$new(
      model = mod_rf,
      data = newdata,  # Original data (no preprocessing applied)
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

# Ensure all grouping columns are atomic vectors
ale_combined <- ale_combined %>%
  mutate(
    predictor = as.character(predictor),  # Convert to character vector
    value = as.numeric(value)             # Convert to numeric vector
  )
ale_dep_nw <- ale_combined

agg_ale <- aggregate(.value ~ predictor + value, 
                     data = ale_combined, 
                     FUN = function(x) c(
                       Mean = mean(x),
                       Lower95 = quantile(x, 0.025),
                       Upper95 = quantile(x, 0.975)
                     ))

# Convert matrix columns to regular columns
agg_ale <- do.call(data.frame, agg_ale)
names(agg_ale)[3:5] <- c("Mean", "Lower95", "Upper95")

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
# Get unique predictors
unique_predictors <- unique(ale_combined$predictor)

# Initialize list to store results
diff_distributions <- list()

# Loop through each predictor
for(pred in unique_predictors) {
  # Subset data for current predictor
  pred_data <- ale_combined[ale_combined$predictor == pred, ]
  
  # Get unique replicates
  replicates <- unique(pred_data$replicate)
  
  # Initialize vector for differences
  diffs <- numeric(length(replicates))
  
  # Calculate differences for each replicate
  for(i in seq_along(replicates)) {
    rep_data <- pred_data[pred_data$replicate == replicates[i], ]
    
    # Get values for level 1 and 2
    value1 <- rep_data$.value[rep_data$value == 1]
    value2 <- rep_data$.value[rep_data$value == 2]
    
    # Store difference
    diffs[i] <- value2 - value1
  }
  
  # Calculate distribution metrics
  diff_distributions[[pred]] <- data.frame(
    predictor = pred,
    median = median(diffs),
    q25 = quantile(diffs, 0.25),
    q75 = quantile(diffs, 0.75),
    lower95 = quantile(diffs, 0.025),
    upper95 = quantile(diffs, 0.975)
  )
}

# Combine all results
agg_diff <- do.call(rbind, diff_distributions)

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


# Convert predictor codes to descriptive labels
agg_diff$predictor <- factor(agg_diff$predictor,
                             levels = names(new_labels),
                             labels = new_labels)

# Order predictors by descending median
agg_diff <- agg_diff[order(-agg_diff$median), ]

# Convert to factor with ordered levels
agg_diff$predictor <- factor(agg_diff$predictor, 
                             levels = agg_diff$predictor)

# Create plot with ordered predictors and median labels
ggplot(agg_diff, aes(x = predictor)) +
  geom_boxplot(
    aes(ymin = lower95, lower = q25, 
        middle = median, upper = q75, ymax = upper95),
    stat = "identity",
    fill = "lightblue",
    width = 0.5
  ) +
  geom_text(
    aes(y = upper95 + 0.005, 
        label = sprintf("%.1f%%", median * 100)),  # Convert to percentage
    hjust = 0,
    size = 5,  # Increased size for median value labels
    color = "darkblue",
    family = "Times New Roman"  # Set font to Times New Roman
  ) +
  coord_flip() +
  labs(x = "Predictor", y = "ALE",
       title = "Effect size of chest pain characterization on CAD prevalence: Depression Cohort") +
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman", size = 14),  # Set font and base size
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Increase title size
    axis.title = element_text(size = 14),  # Increase axis title size
    axis.text = element_text(size = 12),  # Increase axis text size
    panel.grid.major.y = element_blank()
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1)),
    labels = scales::percent_format()  # Convert y-axis to percentages
  )


##################################################
#################################################
#################################################

### Run meta-model for depression cohort on CAD
vars <- c("AGE_BIN", "RIAGENDR", "RIDRETH1", "SMQ020",
          "SMQ040", "HIQ011", "BPQ020", "BMXBMI", "BPQ080",
          "ALQ130", "PAQMV", "DIDTOT","DMDBORNT",
          "INC3", "DMDEDUC2", "SDDSRVYR", "DUQTOT", "CDQ008",
          "CDQ010",
          regions)
outcome <- "CAD"
formula_new <- as.formula(paste(outcome, "~", paste(c(vars), collapse="+")))


all_m0_nw <- list()
tune_grid <- expand.grid(
  mtry = c(8, 16),  # Try different values for mtry
  min.node.size = c(1,5,10),  # Try different values for min.node.size4
  splitrule = "gini"  # Use "gini" for classification (you can also try "extratrees")
)

# Define a custom summary function
customSummary <- function(data, lev = NULL, model = NULL) {
  # 'data$obs' are the observed classes and 'data$pred' are the predicted classes
  sens <- caret::sensitivity(data$pred, data$obs, positive = lev[1])
  spec <- caret::specificity(data$pred, data$obs, negative = lev[2])
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
  tempdata <- all_data[[rep]]
  newdata <- tempdata[tempdata$DEPR_NEW==0,]
  
  # Train a Random Forest model with cross-validation and custom tuning grid
  rf_model <- caret::train(formula_new, 
                           data = newdata, 
                           method = "ranger", 
                           importance = "permutation",
                           replace = TRUE,
                           trControl = cv_control, 
                           num.trees = 100, 
                           tuneGrid = tune_grid,
                           metric = "Specificity"  # Optimize based on your custom 
  )
  all_m0_nw[[rep]] <- rf_model
  
}

# Initialize vectors (one element per model) to store the 0.632 bootstrapped metrics
n_models <- length(all_m0_nw)
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
  model <- all_m0_nw[[i]]
  
  tempdata <- all_data[[i]]
  newdata <- tempdata[tempdata$DEPR_NEW==0,]
  new_dat<-newdata
  levels(new_dat)<-c("0","1")
  oob_data <- newdata[newdata$weights == 0, ]
  ib_data <- newdata[newdata$weights > 0, ]
  
  #### Out-of-Bag (OOB) Performance ####
  # Predicted probabilities and classes on OOB data
  pred_probs_oob <- predict(model, newdata = oob_data, type = "prob")[,2]
  pred_class_oob <- predict(model, newdata = oob_data)
  levels(pred_class_oob)<-c(0,1)
  obs_oob <- oob_data$CAD
  levels(obs_oob)<-c(0,1)
  
  # Determine positive class (assuming the second level)
  positive_class <- levels(obs_oob)[2]
  
  # AUC for OOB
  roc_obj_oob <- pROC::roc(obs_oob, pred_probs_oob)
  auc_oob <- as.numeric(pROC::auc(roc_obj_oob))
  
  # Confusion matrix and related metrics for OOB
  cm_oob <- table(pred_class_oob, obs_oob)
  sens_oob <- cm_oob["1", "1"] / sum(cm_oob["1", ])  # TP / (TP + FN)
  spec_oob <- cm_oob["0", "0"] / sum(cm_oob["0", ])  # TN / (TN + FP)
  
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
  levels(pred_class_ib)<-c(0,1)
  obs_ib <- ib_data$CAD
  levels(obs_ib)<-c(0,1)
  
  # Use the same positive class (assuming levels are consistent)
  positive_class <- levels(obs_ib)[2]
  
  # AUC for IB
  roc_obj_ib <- pROC::roc(obs_ib, pred_probs_ib)
  auc_ib <- as.numeric(pROC::auc(roc_obj_ib))
  
  # Confusion matrix and related metrics for IB
  cm_ib <- table(pred_class_ib, obs_ib)
  sens_ib <- cm_ib["1", "1"] / sum(cm_ib["1", ])  # TP / (TP + FN)
  spec_ib <- cm_ib["0", "0"] / sum(cm_ib["0", ])  # TN / (TN + FP)
  
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


### VIMP

# Create a list to store the importance vectors from each model
importance_list <- lapply(all_m0_nw, function(mod) {
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


#### Sig VIMP
# Filter to keep only variables where the 95% CI does not cross zero
significant_vimp <- importance_summary_sorted[importance_summary_sorted$Lower95 > 0, ]

labels0 <- c("CDQ009F1"="Left Chest Pain",
             "DMDBORNT2"="Not Born in USA",
             "PAQMV1"="No PA",
             "DIDTOT3"="Diabetes: Medication",
             "DUQTOT1"="Use of illicit drugs",
             "ALQ130"="Alcohol use (drinks/day)",
             "SMQ0403"="Current smoker",
             "SMQ0202"="At least 100 cig in lifetime",
             "RIAGENDR2"="Male",
             "DIDTOT4"="Diabetes: Insulin",
             "AGE_BIN5"="> 80 years",
             "AGE_BIN4"="70-79 years",
             "BPQ0202"="History of Hypertension",
             "BPQ0802"="History of Dyslipidemia",
             "CDQ0081"="Severe chest pain for > 30 minutes")


ggplot(significant_vimp, aes(x = reorder(Variable, -Mean), y = Mean)) +
  geom_point(size = 3, color = "red2") +  # Make points larger and color them blue
  geom_errorbar(aes(ymin = Lower95, ymax = Upper95), 
                width = 0.3, 
                size = 1, 
                color = "red2") +  # Thicker error bars in dark red
  coord_flip() +
  theme_minimal(base_family = "Times New Roman", base_size = 16) +  # Set Times New Roman globally
  labs(
    title = "Significant Variable Importance with 95% CI \n for Predicting CAD in No Depression Model",
    x = "Variable", 
    y = "Permuted Performance Difference"
  ) +
  scale_x_discrete(labels = labels0) +  # Ensure labels are properly displayed
  theme(
    axis.text.x = element_text(size = 16, family = "Times New Roman"),  
    axis.text.y = element_text(size = 16, family = "Times New Roman"),  
    axis.title.x = element_text(size = 18, face = "bold", family = "Times New Roman"),  
    axis.title.y = element_text(size = 18, face = "bold", family = "Times New Roman"),  
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5, family = "Times New Roman")  
  )






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
  ale_list[[var]] <- lapply(seq_along(all_m0_nw), function(i) {
    # Extract the ranger model from caret
    mod_rf <- all_m0_nw[[i]]
    tempdata <- all_data[[i]]
    newdata <- tempdata[tempdata$DEPR_NEW==0,]
    
    # Create predictor object with corrected model and data
    predictor_obj <- Predictor$new(
      model = mod_rf,
      data = newdata,  # Original data (no preprocessing applied)
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


# Ensure all grouping columns are atomic vectors
ale_combined <- ale_combined %>%
  mutate(
    predictor = as.character(predictor),  # Convert to character vector
    value = as.numeric(value)             # Convert to numeric vector
  )

ale_nodep_nw <- ale_combined

agg_ale <- aggregate(.value ~ predictor + value, 
                     data = ale_combined, 
                     FUN = function(x) c(
                       Mean = mean(x),
                       Lower95 = quantile(x, 0.025),
                       Upper95 = quantile(x, 0.975)
                     ))

# Convert matrix columns to regular columns
agg_ale <- do.call(data.frame, agg_ale)
names(agg_ale)[3:5] <- c("Mean", "Lower95", "Upper95")

# Plot
ggplot(agg_ale, aes(x = factor(value), y = Mean)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Lower95, ymax = Upper95), width = 0.1) +
  facet_wrap(~ predictor, scales = "free_y") +
  labs(x = "Level", y = "Accumulated Local Effect") +
  theme_minimal()

# For each predictor and replicate, compute the difference between level "1" and level "0".
# (Assuming the levels are coded as "0" and "1"; adjust if necessary.)
# Get unique predictors
unique_predictors <- unique(ale_combined$predictor)

# Initialize list to store results
diff_distributions <- list()

# Loop through each predictor
for(pred in unique_predictors) {
  # Subset data for current predictor
  pred_data <- ale_combined[ale_combined$predictor == pred, ]
  
  # Get unique replicates
  replicates <- unique(pred_data$replicate)
  
  # Initialize vector for differences
  diffs <- numeric(length(replicates))
  
  # Calculate differences for each replicate
  for(i in seq_along(replicates)) {
    rep_data <- pred_data[pred_data$replicate == replicates[i], ]
    
    # Get values for level 1 and 2
    value1 <- rep_data$.value[rep_data$value == 1]
    value2 <- rep_data$.value[rep_data$value == 2]
    
    # Store difference
    diffs[i] <- value2 - value1
  }
  
  # Calculate distribution metrics
  diff_distributions[[pred]] <- data.frame(
    predictor = pred,
    median = median(diffs),
    q25 = quantile(diffs, 0.25),
    q75 = quantile(diffs, 0.75),
    lower95 = quantile(diffs, 0.025),
    upper95 = quantile(diffs, 0.975)
  )
}

# Combine all results
agg_diff <- do.call(rbind, diff_distributions)

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


# Convert predictor codes to descriptive labels
agg_diff$predictor <- factor(agg_diff$predictor,
                             levels = names(new_labels),
                             labels = new_labels)

# Order predictors by descending median
agg_diff <- agg_diff[order(-agg_diff$median), ]

# Convert to factor with ordered levels
agg_diff$predictor <- factor(agg_diff$predictor, 
                             levels = agg_diff$predictor)

# Create plot with ordered predictors and median labels
ggplot(agg_diff, aes(x = predictor)) +
  geom_boxplot(
    aes(ymin = lower95, lower = q25, 
        middle = median, upper = q75, ymax = upper95),
    stat = "identity",
    fill = "lightblue",
    width = 0.5
  ) +
  geom_text(
    aes(y = upper95 + 0.005, 
        label = sprintf("%.1f%%", median * 100)),  # Convert to percentage
    hjust = 0,
    size = 5,  # Increased size for median value labels
    color = "darkblue",
    family = "Times New Roman"  # Set font to Times New Roman
  ) +
  coord_flip() +
  labs(x = "Predictor", y = "ALE",
       title = "Effect size of chest pain characterization on CAD prevalence: No Depression Cohort") +
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman", size = 14),  # Set font and base size
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Increase title size
    axis.title = element_text(size = 14),  # Increase axis title size
    axis.text = element_text(size = 12),  # Increase axis text size
    panel.grid.major.y = element_blank()
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1)),
    labels = scales::percent_format()  # Convert y-axis to percentages
  )
