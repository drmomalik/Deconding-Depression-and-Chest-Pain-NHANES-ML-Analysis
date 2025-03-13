
library(MASS) # For multivariate probit model (example purpose)

### How to run non-survey package analysis with complex survey design
## 1. Create a svydesign object using the imputed version of the original data
## 2. Create a replicate weight using bootstrap and n value of replications
## 3. Subset the replicate weight design matrix
## 4. Create weights and stabilize 
## 5. Find estimate for base data (using pweights)
## 6. Find estimate for n copy of replicates (using repweights)
## 7. Use distribution of estimates to get a bootstrap estimate of variance or CI using percentiles





#Create a svydesign object using the imputed version of the original data

imputation <- complete(imputations, 1)

design <- svydesign(
  id = ~SDMVPSU, 
  weights = ~MEC15YR, 
  strata = ~SDMVSTRA, 
  data = imputation,
  nest = TRUE)

set.seed(123)  
rep_design <- as.svrepdesign(design, type = "bootstrap", replicates=100)
rep_design <- subset(rep_design, !apply(
  rep_design$variables[, c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")],
  1, 
  function(x) all(x == 0 | is.na(x))  # Keep rows where not all values of CDQ009A-CDQ009H are 0 or NA
))

# Create weights and stabilize 
repweights <- weights(rep_design, type = "analysis")
repweights <- repweights/mean(repweights)
pweights <- weights(rep_design, type = "sampling")
pweights <- pweights/mean(pweights)

### Example use of bootstrap replicated weights to determine mean and SE
# Base mean estimate
base_mean <- sum(pweights * rep_design$variables$BMXBMI) / sum(pweights)
# Replicate mean estimates
rep_means <- apply(repweights, 2, function(w) {
  sum(w * rep_design$variables$BMXBMI) / sum(w)
})
# Variance of the mean estimate
rep_var <- mean((rep_means - base_mean)^2)

# Standard error
rep_se <- sqrt(rep_var)

# Print results
cat("Base Mean:", base_mean, "\n")
cat("Standard Error (Replicate Variance):", rep_se, "\n")


### Example use of boostrap replicated weights to determine simple glm 

# Base model using pweights
base_model <- glm(CDQ009A ~ -1+DEPR_BIN + AGE_BIN + RIAGENDR + RIDRETH1 + SMQ040 + 
                    HIQ011 + BPQ020 + BMXBMI + BPQ080 + ALQ130 + MEDDEP + DMDBORNT + 
                    PAQMV + CADTOT + DIDTOT + DUQTOT + INC_BIN, family = binomial, data = rep_design$variables, weights = pweights)
base_coef <- coef(base_model)

# Replicate estimates (use apply for each replicate weight column)
rep_coefs <- apply(repweights, 2, function(w) {
  model <- glm(CDQ009A ~ -1+DEPR_BIN + AGE_BIN + RIAGENDR + RIDRETH1 + SMQ040 + 
                 HIQ011 + BPQ020 + BMXBMI + BPQ080 + ALQ130 + MEDDEP + DMDBORNT + 
                 PAQMV + CADTOT + DIDTOT + DUQTOT + INC_BIN, family = binomial, data = rep_design$variables, weights = w)
  coef(model)  # Extract coefficients
})

# Transpose rep_coefs to match coefficient names and make it easier to compute variance
rep_coefs <- t(rep_coefs)

# Variance for each coefficient (compare each replicate coefficient with the base model coefficient)
rep_var <- apply(rep_coefs, 2, function(coef) {
  mean((coef - base_coef)^2)
})

# Standard errors
rep_se <- sqrt(rep_var)

# Combine results into a data frame with 5 coefficients and 5 standard errors
coef_summary <- data.frame(
  Coefficient = base_coef,
  StdError = rep_se
)

# Print results
print(coef_summary)



### Penalized Log Reg
# Load glmnet library
library(glmnet)
library(pROC)
library(survey)

imputation <- complete(imputations, 1)

design <- svydesign(
  id = ~SDMVPSU, 
  weights = ~MEC15YR, 
  strata = ~SDMVSTRA, 
  data = imputation,
  nest = TRUE)

set.seed(123)  
set_design <- as.svrepdesign(design, type = "bootstrap", replicates=100)
rep_design <- subset(set_design, !apply(
  set_design$variables[, c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")],
  1, 
  function(x) all(x == 0 | is.na(x))  # Keep rows where not all values of CDQ009A-CDQ009H are 0 or NA
))




# Create weights and stabilize 
repweights <- weights(rep_design, type = "analysis")
repweights <- repweights/mean(repweights)
pweights <- weights(rep_design, type = "sampling")
pweights <- pweights/mean(pweights)


data <- rep_design$variables
data$pweights <- pweights
data$CDQ81 <- if_else(data$CDQ008 == 1 | data$CDQ010 == 1, 1, 0)


# Prepare the data for glmnet
# Convert predictor variables to matrix form
X <- model.matrix(~~ AGE_BIN + RIAGENDR + RIDRETH1 +
                    HIQ011 + I(BMXBMI^2) + INC3 + DMDEDUC2 + 
                    PAQMV + I(ALQ130^2) + DIDTOT + DUQTOT + SMQ040
                  + DMDBORNT + BPQ020 + BPQ080 + SDDSRVYR + CDQ008 +
                    DEPR_BIN*MEDDEP + DEPR_BIN*CADTOT + CDQ010,
                  data = rep_design$variables)


# Convert the response variable to a vector
y <- rep_design$variables$CDQ009B

# Create penalty factors: lower penalty for priority_vars
priority_vars <- c("DEPR_BIN2", "DEPR_BIN2:MEDDEP1", "DEPR_BIN2:CADTOT1")
penalty_factors <- ifelse(colnames(X) %in% priority_vars, 0, 1)


# Define the base model using glmnet
base_model <- cv.glmnet(X, y, family = "binomial", weights = pweights, 
                        penalty.factor = penalty_factors,
                         alpha = 0, nfolds = 10, type.measure = "auc"
                        , relax=TRUE )



coef(base_model, s="lambda.min")

#Assess ROC-AUC
lambda_min <- base_model$lambda.min
auc <- base_model$cvm[base_model$lambda == lambda_min]
print(auc)
plot(base_model)



###
mod <- glm(CDQ009H ~ AGE_BIN + RIAGENDR + RIDRETH1 +
             HIQ011 + I(BMXBMI^2) + INC3 + DMDEDUC2 + 
             PAQMV + I(ALQ130^2) + DIDTOT + DUQTOT + SMQ040*SMQ020
           + DMDBORNT + BPQ020 + BPQ080 + SDDSRVYR +
             DEPR_LVL*MEDDEP + DEPR_LVL*CADTOT + 
           CDQ008 + CDQ010,
           data = rep_design$variables, weights=pweights, family=quasibinomial(link="logit"))

pred_probs <- predict(mod, type = "response")

y <- rep_design$variables$CDQ009H

# Compute ROC curve
roc_curve <- roc(y, pred_probs)

# Calculate AUC
auc_value <- auc(roc_curve)
print(auc_value)

# Optionally, plot ROC curve
plot(roc_curve)



# Install the superml package if needed
# install.packages("superml")

library(superml)

data <- rep_design$variables
data$pweights <- pweights

# Prepare your formula and data
formula <- CDQ009B ~ AGE_BIN + RIAGENDR + RIDRETH1 + HIQ011 + BMXBMI + INC_BIN + DMDEDUC2 +
  PAQMV + ALQ130 + DIDTOT + DUQTOT + SMQ040 + DMDBORNT + BPQ020 + BPQ080 + SDDSRVYR +
  DEPR_LVL + MEDDEP + CADTOT

# Fit logistic regression with cross-validation and weights
cv_model_superml <- train(
  formula = formula,
  data = data,
  model_type = "logistic",  # Specify logistic regression
  kfold = 10,  # Number of cross-validation folds
  weights = pweights,  # Include weights
  metric = "auc"  # Optimize for AUC
)

# View results
cv_model_superml$results

###



## Find marginal effect of DEPR_BIN 
probs <- table(rep_design$variables$MEDDEP, rep_design$variables$CADTOT) / nrow(rep_design$variables)
marginal_effect <- coef(base_model, s="lambda.min")@x[2] + coef(base_model, s="lambda.min")@x[48]*probs["1","0"]
+ coef(base_model, s="lambda.min")@x[49]*probs["0","1"] + 
  (coef(base_model, s="lambda.min")@x[48]+coef(base_model, s="lambda.min")@x[49])*probs["1","1"]
print(marginal_effect)

plot(base_model)
mar = c(4.5, 4.5, 4, 1)

# Extract coefficients for the base model and convert them to a dense vector
base_coef <- coef(base_model)  # Use optimal lambda
base_coef_dense <- as.vector(base_coef)

# Replicate estimates (use apply for each replicate weight column)
rep_coefs <- apply(repweights, 2, function(w) {
  # Fit the penalized logistic regression model for each replicate weight
  model <- glmnet(X, y, family = "binomial", weights = w, alpha = 1)
  coef(model)  # Extract coefficients for the replicate model
})

# Check if rep_coefs is a list and convert it to a matrix if so
if (is.list(rep_coefs)) {
  rep_coefs <- do.call(cbind, rep_coefs)  # Convert list to matrix by binding columns
}

# Ensure rep_coefs is a numeric matrix (dense format)
rep_coefs_dense <- as.matrix(rep_coefs)  # Convert to dense matrix if necessary

# Variance for each coefficient (compare each replicate coefficient with the base model coefficient)
rep_var <- apply(rep_coefs_dense, 1, function(coef) {
  mean((coef - base_coef_dense)^2)
})

# Standard errors
rep_se <- sqrt(rep_var)

# Combine results into a data frame with coefficients and standard errors
coef_summary <- data.frame(
  Coefficient = base_coef_dense,
  StdError = rep_se
)

# Print results
print(coef_summary)


#### 

# Load necessary libraries
library(lme4)
library(dplyr)
library(tidyr)

imputation <- complete(imputations, 1)

# Reshape the data to long format
long_data <- imputation %>%
  pivot_longer(cols = starts_with("CDQ009"), 
               names_to = "outcome_name", 
               values_to = "outcome") %>%
  mutate(outcome_name = as.factor(outcome_name)) 

design <- svydesign(
  id = ~SDMVPSU, 
  weights = ~MEC15YR, 
  strata = ~SDMVSTRA, 
  data = long_data,
  nest = TRUE)

set.seed(123)  
rep_design <- as.svrepdesign(design, type = "bootstrap", replicates=100)
# Remove SEQNs where all outcomes are 0

# Example: Assuming `rep_design` is your svyrep.design object
# Step 1: Extract the data and replicate weights
rep_data <- rep_design$variables
rep_weights <- as.matrix(rep_design$repweights)
pweights <- as.vector(rep_design$pweights)

# Step 2: Identify SEQNs to keep based on the filtering logic
valid_seqn <- rep_data %>%
  group_by(SEQN) %>%
  filter(any(outcome != 0)) %>%
  pull(SEQN) %>%
  unique()

# Step 3: Filter the dataset and replicate weights
filtered_rep_data <- rep_data %>%
  filter(SEQN %in% valid_seqn)

filtered_rep_weights <- rep_weights[rep_data$SEQN %in% valid_seqn, ]

filtered_pweights <- pweights[rep_data$SEQN %in% valid_seqn]

# Step 4: Recreate the svyrep.design object with the filtered data
filtered_rep_design <- svrepdesign(
  variables = filtered_rep_data,
  weights = filtered_pweights, # Match weights to filtered data
  repweights = filtered_rep_weights,                            # Match replicate weights
  scale = rep_design$scale,
  rscales = rep_design$rscales,
  type = rep_design$type
)

# The filtered survey design is now ready




# Create weights and stabilize 
repweights <- weights(filtered_rep_design, type = "analysis")
repweights <- repweights/mean(repweights)
pweights <- weights(filtered_rep_design, type = "sampling")
pweights <- pweights/mean(pweights)



# Fit the GLMM with random intercepts for SEQN
model <- glmer(
  outcome ~ DEPR_BIN + outcome_name + (1 | SEQN),  # Fixed effect for DEPR_BIN and outcome type
  family = binomial(),                             # Binary response
  weights = pweights,                              # Incorporate weights
  data = filtered_rep_design$variables
)

library(brms)
model <- brm(
  outcome ~ DEPR_LVL + outcome_name + (1 | SEQN),
  family = bernoulli(),
  weights
  data = filtered_rep_design$variables,
  chains = 4, iter = 2000
)




###glmmTMB

# Load glmnet library
library(glmnet)
library(pROC)
library(survey)
library(glmmTMB)


imputation <- complete(imputations, 1)


design <- svydesign(
  id = ~SDMVPSU, 
  weights = ~MEC15YR, 
  strata = ~SDMVSTRA, 
  data = imputation,
  nest = TRUE)

set.seed(123)  
set_design <- as.svrepdesign(design, type = "bootstrap", replicates=100)
rep_design <- subset(set_design, !apply(
  set_design$variables[, c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")],
  1, 
  function(x) all(x == 0 | is.na(x))  # Keep rows where not all values of CDQ009A-CDQ009H are 0 or NA
))



# Create weights and stabilize 
repweights <- weights(rep_design, type = "analysis")
repweights <- repweights/mean(repweights)
pweights <- weights(rep_design, type = "sampling")
pweights <- pweights/mean(pweights)


# Extract data frame and add back pweights
data <- rep_design$variables
data$pweights <- pweights

# Reshape the data to long format
long_data <- data %>%
  pivot_longer(cols = starts_with("CDQ009"), # Select columns CDQ009A to CDQ009H
               names_to = "outcome_id",      # New column to hold response IDs
               values_to = "response") %>%   # Column to hold response values
  mutate(outcome_id = factor(outcome_id))    # Convert outcome_id to a factor


# Define the model formula
formula <- response ~ DEPR_BIN*MEDDEP + DEPR_BIN*CADTOT +
  (1 | SEQN) + (1 | outcome_id)  # Random effects for subjects and outcomes

# Fit the model
model <- glmmTMB(
  formula,
  family = binomial(link = "logit"),  # Logistic regression
  weights = pweights,                # Observation weights
  data = long_data                  # Use reshaped long-format data
)

# Summary of the model
summary(model)

# Extract random effects (intercepts) for individuals (SEQN) and outcomes (outcome_id)
random_effects <- ranef(model)

# View random intercepts for individuals
individual_random_effects <- random_effects$cond$SEQN

# View random intercepts for outcomes
outcome_random_effects <- random_effects$cond$outcome_id

# Combine the random effects into a single data frame for analysis
random_effects_df <- data.frame(individual_random_effects, outcome_random_effects)

# Compute the correlation matrix between the random intercepts of outcomes
cor_matrix <- cor(outcome_random_effects)

# Print the correlation matrix between the outcomes
print(cor_matrix)



###



