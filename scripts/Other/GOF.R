
### GOF 

# AIC
# Extract log-likelihoods and parameters
logLik_mod <- logLik(mod)
logLik_mod2 <- logLik(mod2)

k_mod <- length(coef(mod))
k_mod2 <- length(coef(mod2))

# Compute AIC
aic_mod <- -2 * logLik_mod + 2 * k_mod
aic_mod2 <- -2 * logLik_mod2 + 2 * k_mod2

c(AIC_mod = aic_mod, AIC_mod2 = aic_mod2)
## 

plot(base_model)
mar = c(4.5, 4.5, 4, 1)

coef(base_model, s="lambda.min")

#Assess ROC-AUC
# Predicted probabilities for the optimal lambda (lambda.min)
pred_probs <- predict(mod, type = "response")

y <- rep_design$variables$CDQ009B
  
# Compute ROC curve
roc_curve <- roc(y, pred_probs)

# Calculate AUC
auc_value <- auc(roc_curve)
print(auc_value)

# Optionally, plot ROC curve
plot(roc_curve)


# Deviance residuals
deviance_residuals <- sqrt(-2 * (as.numeric(y) * log(pred_probs) + (1 - as.numeric(y)) * log(1 - pred_probs)))

# Pearson residuals
pearson_residuals <- (as.numeric(y) - pred_probs) / sqrt(pred_probs * (1 - pred_probs))

# Identify indexes with deviance residuals >= 2
outlier_indexes <- which(deviance_residuals >= 2)

# View the indexes
print(outlier_indexes)


### AUC for glmmTMB

library(pROC)

# Extract predicted probabilities from the glmmTMB model
pred_probs <- predict(model, type = "response")  # Predicted probabilities

# Initialize a data frame to store AUC results for each outcome
AUC_results <- data.frame(outcome_id = unique(long_data$outcome_id), AUC = NA)

# Loop through each outcome and compute AUC
for (outcome in unique(long_data$outcome_id)) {
  # Get the subset of data for the current outcome
  outcome_data <- subset(long_data, outcome_id == outcome)
  
  # Extract the predicted probabilities for the current outcome
  # Ensure the order of predicted probabilities matches the subset of data
  outcome_pred_probs <- pred_probs[long_data$outcome_id == outcome]
  
  # Compute AUC for each outcome
  roc_curve <- roc(outcome_data$response, outcome_pred_probs)
  AUC_results$AUC[AUC_results$outcome_id == outcome] <- auc(roc_curve)
}

# Print AUC results for each outcome
print(AUC_results)

###

### Dev Res for glmmTMB

# Initialize a data frame to store deviance residuals for each outcome
deviance_results <- data.frame(matrix(ncol = length(unique(long_data$outcome_id)), 
                                      nrow = nrow(long_data)))

# Name the columns as the unique outcomes
colnames(deviance_results) <- unique(long_data$outcome_id)

# Initialize a list to store indices for outliers (deviance residuals >= 2)
outlier_indices_list <- list()

# Loop through each outcome and compute deviance residuals
for (outcome in unique(long_data$outcome_id)) {
  # Get the subset of data for the current outcome
  outcome_data <- subset(long_data, outcome_id == outcome)
  
  # Extract the predicted probabilities for the current outcome
  outcome_pred_probs <- pred_probs[long_data$outcome_id == outcome]
  
  # Compute deviance residuals for all observations
  deviance_residuals <- residuals(model, type = "deviance")
  
  # Store the deviance residuals in the correct column for the current outcome
  deviance_results[[outcome]] <- deviance_residuals[long_data$outcome_id == outcome]
  
  # Plot deviance residuals vs fitted values for the current outcome
  plot <- ggplot(data.frame(fitted = outcome_pred_probs, 
                            deviance_residuals = deviance_residuals[long_data$outcome_id == outcome]), 
                 aes(x = fitted, y = deviance_residuals)) +
    geom_point(color = "blue", alpha = 0.5) +
    geom_smooth(method = "loess", color = "red", se = FALSE) +
    theme_minimal() +
    labs(title = paste("Deviance Residuals vs Fitted Values for Outcome", outcome),
         x = "Fitted Values (Predicted Probabilities)", y = "Deviance Residuals")
  
  # Explicitly print the plot for the current outcome
  print(plot)
  
  # Identify indices where deviance residuals are greater than or equal to 2
  outlier_indices <- which(deviance_residuals >= 2 & long_data$outcome_id == outcome)
  
  # Store the indices in the outlier_indices_list
  outlier_indices_list[[outcome]] <- outlier_indices
}

# Print deviance results for each outcome
print(deviance_results)

# Print the outlier indices for each outcome
print(outlier_indices_list)


###





# Plot residuals vs fitted values
plot(pred_probs, deviance_residuals, 
     main = "Deviance Residuals vs Fitted Values", 
     xlab = "Fitted Values", 
     ylab = "Deviance Residuals")
abline(h = 0, col = "red")

hist(deviance_residuals, main = "Histogram of Deviance Residuals", 
     xlab = "Deviance Residuals", breaks = 30)

# Q-Q plot
qqnorm(deviance_residuals, main = "Q-Q plot of Deviance Residuals")
qqline(deviance_residuals, col = "red")

library(car)
vif(base_model$glmnet.fit)



### Defining outliers for each model (CDQ009A:H)


# Response variables
responses <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")

# Data frame to store outlier indexes
outliers <- data.frame(response = responses, outlier_indexes = I(vector("list", length(responses))))

for (response in responses) {
  # Prepare the data for glmnet
  # Convert predictor variables to matrix form
  X <- model.matrix(~ RIAGENDR*AGE_BIN + RIDRETH1 + SMQ020 +
                      HIQ011 + BPQ020 + BPQ080 +  DMDBORNT + 
                       DIDTOT + INC_BIN + DEPR_BIN*MEDDEP + DEPR_BIN*CADTOT
                    + SMQ040*ALQ130*DUQTOT + PAQMV*BMI_LVL,
                    data = rep_design$variables)
  
  # Convert the response variable
  y <- rep_design$variables[[response]]
  
  # Create penalty factors: lower penalty for priority_vars
  priority_vars <- c("DEPR_BIN2", "DEPR_BIN2:MEDDEP1", "DEPR_BIN2:CADTOT1")
  penalty_factors <- ifelse(colnames(X) %in% priority_vars, 0, 1)
  
  # Define the base model using glmnet
  base_model <- cv.glmnet(X, y, family = "binomial", weights = rep_design$pweights,
                          penalty.factor = penalty_factors, alpha = 0, nfolds = 10, 
                          type.measure = "auc", standardize = TRUE)
  
  # Assess ROC-AUC
  # Predicted probabilities for the optimal lambda (lambda.min)
  pred_probs <- predict(base_model, newx = X, s = "lambda.min", type = "response")
  
  # Convert y to numeric
  y <- as.numeric(y)
  
  # Deviance residuals
  deviance_residuals <- sqrt(-2 * (y * log(pred_probs) + (1 - y) * log(1 - pred_probs)))
  
  # Identify indexes with deviance residuals >= 2
  outlier_indexes <- which(deviance_residuals >= 3)
  
  # Append outlier indexes to the outliers data frame
  outliers$outlier_indexes[outliers$response == response] <- list(outlier_indexes)
}

# Print the resulting outliers data frame
print(outliers)





