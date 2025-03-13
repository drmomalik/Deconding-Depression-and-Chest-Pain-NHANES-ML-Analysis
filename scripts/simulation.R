# ----------------------------
# 1. Load Required Packages
# ----------------------------
library(survey)         # For survey design and replicate weights
library(caret)          # For model training and cross validation
library(randomForest)   # For the RF algorithm (used by caret via method = "rf")
library(iml)            # For computing Accumulated Local Effects (ALE)
library(dplyr)          # For data manipulation
library(ggplot2)        # For plotting
library(purrr)          # For list manipulation
set.seed(123)

# ----------------------------
# 2. Simulate Population Data with a Known Treatment Effect
# ----------------------------
# We create a population where a binary predictor "trt" increases the probability of a binary outcome
# by approximately 15 percentage points at a typical value of a continuous covariate (X).
# For example, we set up a logistic model so that when X = 0:
#   P(Y=1 | trt=0) ≈ 0.20  and  P(Y=1 | trt=1) ≈ 0.35.


### STEP 1: Simulate the Superpopulation

N_pop <- 10000

# Simulate covariates
pop_trt <- rbinom(N_pop, 1, 0.4)     # Treatment indicator
pop_X   <- rnorm(N_pop)              # Continuous covariate X
pop_Z   <- rnorm(N_pop)              # Additional covariate Z for stratification

# Set true parameters
beta0 <- -1.3863        # Intercept (baseline logit corresponding roughly to 20%)
beta1 <- 0.7673         # Treatment effect (~15% absolute change at reference)
beta2 <- 0.3            # Effect of X
beta3 <- 0.5            # Effect of Z
sigma_cluster <- 0.2    # Cluster-level variability

# Create clusters and strata:
# Form strata based on quantiles of Z
strata_levels <- paste0("stratum_", 1:100)
pop_strata <- cut(pop_Z, breaks = quantile(pop_Z, probs = seq(0, 1, length.out = 101)),
                  include.lowest = TRUE, labels = strata_levels)

# Within each stratum, assign clusters randomly from 1 to 10
pop_cluster <- ave(as.numeric(pop_strata), pop_strata,
                   FUN = function(x) paste0("cluster_", sample(1:10, length(x), replace = TRUE)))

# Create cluster-level random effects
unique_clusters <- unique(pop_cluster)
cluster_effects <- rnorm(length(unique_clusters), mean = 0, sd = sigma_cluster)
names(cluster_effects) <- unique_clusters

# Compute the linear predictor (including cluster effect)
linpred <- beta0 + beta1 * pop_trt + beta2 * pop_X + beta3 * pop_Z +
  cluster_effects[pop_cluster]
pop_prob <- exp(linpred) / (1 + exp(linpred))
pop_Y <- rbinom(N_pop, 1, pop_prob)

pop_data <- data.frame(
  Y = pop_Y,
  trt = pop_trt,
  X = pop_X,
  Z = pop_Z,
  cluster = pop_cluster,
  strata = pop_strata
)

### STEP 2: Impose a Complex Survey Sampling Design

# Define selection probabilities that depend on treatment (and possibly X)
pop_data <- pop_data %>%
  mutate(p_select = ifelse(trt == 1, 0.8, 0.5),
         # Optionally modify selection further by X to add extra imbalance:
         p_select = p_select * (1 + 0.2 * (X - mean(X))),
         # Keep probabilities within a reasonable range:
         p_select = pmin(pmax(p_select, 0.1), 1))

# Draw a survey sample based on these probabilities
pop_data$selected <- rbinom(N_pop, 1, pop_data$p_select)
survey_data <- pop_data %>% filter(selected == 1)
survey_data <- survey_data %>% mutate(weight = 1 / p_select)

# For modeling with caret we need factors.
survey_data$Y_factor    <- factor(survey_data$Y, levels = c(0, 1), labels = c("No", "Yes"))
survey_data$trt_factor  <- factor(survey_data$trt, levels = c(0, 1), labels = c("Control", "Treatment"))


# ----------------------------
# 4. Create Survey Design and Replicate Weights
# ----------------------------
design_obj <- svydesign(ids = ~cluster,
                        strata = ~strata,
                        weights = ~weight,
                        data = survey_data,
                        nest = TRUE)

# Create replicate weights via bootstrap; here we use 50 replicates.
rep_design <- as.svrepdesign(design_obj, type = "bootstrap", replicates = 100)
repweights <- weights(rep_design, type = "analysis")
repweights <- repweights/mean(repweights) # Stabilize case weights 

# Each column of repweights is a replicate (accessible as repweights[, i]).

# ----------------------------
# 5. Set Up caret Training Control and Model Formula
# ----------------------------
tune_grid <- expand.grid(
  mtry = c(1, 2),  # Try different values for mtry
  min.node.size = c(1, 5, 10),  # Try different values for min.node.size
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


# Our model will predict the binary outcome from treatment and the continuous predictor.
formula_rf <- Y_factor ~ trt_factor + X + Z

# ----------------------------
# 6. Define a Helper Function to Compute ALE Difference for the Treatment Effect
# ----------------------------
# This function fits an iml Predictor on a given model, computes ALE for "trt_factor",
# and returns the difference: ALE(Treatment) minus ALE(Control).
compute_ALE_diff <- function(model, data) {
  # Custom prediction function for caret rf model: returns probability for "Yes"
  pred_fun <- function(model, newdata) {
    predict(model, newdata = newdata, type = "prob")[, "Yes"]
  }
  predictor_obj <- Predictor$new(
    model = model,
    data = survey_data[, c("trt_factor", "X", "Z")],
    y = survey_data$Y_factor,
    predict.function = pred_fun
  )
  # Compute ALE for "trt_factor" with grid.size = 2 (since it's binary)
  ale_obj <- FeatureEffect$new(
    predictor = predictor_obj,
    feature = "trt_factor",
    method = "ale",
    grid.size = 2
  )
  ale_data <- ale_obj$results
  # Here we assume that the ALE results include rows labeled by trt_factor (e.g., "Control" and "Treatment").
  effect_diff <- ale_data$.value[ale_data$trt_factor == "Treatment"] -
    ale_data$.value[ale_data$trt_factor == "Control"]
  return(effect_diff)
}

# ----------------------------
# 7. Fit Weighted Random Forest Models Using Replicate Weights and Compute ALE Differences
# ----------------------------
B <- ncol(repweights)   # Number of replicates (50)
all_model_weighted <- list()      # To store models for each replicate
ale_diff_weighted <- numeric(B)   # To store ALE difference from each replicate


for (i in 1:B) {
  cat("Fitting weighted model for replicate", i, "\n")
  set.seed(100 + i)
  model_i <- caret::train(formula_rf, 
                   data = survey_data, 
                   method = "ranger", 
                   weights = repweights[,i],
                   trControl = cv_control,
                   replace = TRUE, 
                   num.trees = 50,         
                   tuneGrid = tune_grid,
                   metric= "auc")
  all_model_weighted[[i]] <- model_i
  ale_diff_weighted[i] <- compute_ALE_diff(model_i, survey_data)
  cat("  Weighted replicate", i, "ALE diff:", round(ale_diff_weighted[i], 3), "\n")
}

# ----------------------------
# 8. Fit Unweighted Random Forest Models (Ignoring Weights) and Compute ALE Differences
# ----------------------------
all_model_unweighted <- list()    # To store unweighted models
ale_diff_unweighted <- numeric(B) # To store ALE differences for each replicate

for (i in 1:B) {
  cat("Fitting unweighted model for replicate", i, "\n")
  set.seed(200 + i)
  model_i <- caret::train(formula_rf, 
                   data = survey_data, 
                   method = "rf", 
                   importance = TRUE,
                   weights = rep(1, nrow(survey_data)),  # Ignore weights: use constant weights
                   trControl = cv_control,
                   ntree = 50,
                   nodesize = 10,
                   tuneLength = 3)
  all_model_unweighted[[i]] <- model_i
  ale_diff_unweighted[i] <- compute_ALE_diff(model_i, survey_data)
  cat("  Unweighted replicate", i, "ALE diff:", round(ale_diff_unweighted[i], 3), "\n")
}

# ----------------------------
# 9. Compare the ALE Difference Distributions
# ----------------------------
# Combine the results into a data frame for easy comparison.
ale_diff_df <- data.frame(
  Method   = rep(c("Weighted", "Unweighted"), each = B),
  ALE_diff = c(ale_diff_weighted, ale_diff_unweighted)
)
print(ale_diff_df)

# Plot histograms of the ALE differences.
ggplot(ale_diff_df, aes(x = ALE_diff, fill = Method)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 15) +
  geom_vline(xintercept = 0.15, linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = 0.15, y = max(table(cut(ale_diff_df$ALE_diff, 15))), 
           label = "True Population Effect", color = "red", hjust = -0.1, vjust = 1) +
  labs(title = "Distribution of ALE differences using\nsimulated complex survey data (Treatment - Control)",
       x = "ALE Difference", y = "Count") +
  theme_minimal()



# Summarize the ALE differences.
weighted_summary <- data.frame(
  Method   = "Weighted",
  Mean     = mean(ale_diff_weighted),
  SD       = sd(ale_diff_weighted),
  Lower95  = quantile(ale_diff_weighted, 0.025),
  Upper95  = quantile(ale_diff_weighted, 0.975)
)
unweighted_summary <- data.frame(
  Method   = "Unweighted",
  Mean     = mean(ale_diff_unweighted),
  SD       = sd(ale_diff_unweighted),
  Lower95  = quantile(ale_diff_unweighted, 0.025),
  Upper95  = quantile(ale_diff_unweighted, 0.975)
)
summary_all <- rbind(weighted_summary, unweighted_summary)
print(summary_all)

###########################
##########################

# Run sim multiple times to get 95% coverage probability 

# Load required libraries
library(caret)
library(ranger)
library(randomForest)
library(iml)
library(ggplot2)
library(survey)

#############################
# PRE-DEFINED CODE SETTINGS
#############################

# 1. Define your tune grid, custom summary function, and CV control
tune_grid <- expand.grid(
  mtry = c(1, 2),
  min.node.size = c(1, 5, 10),
  splitrule = "gini"
)

customSummary <- function(data, lev = NULL, model = NULL) {
  sens <- sensitivity(data$pred, data$obs, positive = lev[1])
  spec <- specificity(data$pred, data$obs, negative = lev[2])
  gmean <- sqrt(sens * spec)
  out <- c(Sensitivity = sens, Specificity = spec, GMean = gmean)
  return(out)
}

cv_control <- trainControl(method = "cv", 
                           number = 5, 
                           classProbs = TRUE,
                           summaryFunction = customSummary)

# 2. Define the model formula (assumes your survey_data has these columns)
formula_rf <- Y_factor ~ trt_factor + X + Z

# 3. Define the helper function to compute the ALE difference
compute_ALE_diff <- function(model, data) {
  # Custom prediction function: returns probability for "Yes"
  pred_fun <- function(model, newdata) {
    predict(model, newdata = newdata, type = "prob")[, "Yes"]
  }
  
  # Here we use all the predictors that the model was trained on.
  predictor_obj <- Predictor$new(
    model = model,
    data = data[, c("trt_factor", "X", "Z")],
    y = data$Y_factor,
    predict.function = pred_fun
  )
  
  # Compute ALE for "trt_factor" with grid.size = 2 (for binary feature)
  ale_obj <- FeatureEffect$new(
    predictor = predictor_obj,
    feature = "trt_factor",
    method = "ale",
    grid.size = 2
  )
  ale_data <- ale_obj$results
  
  # Assume that the ALE results include rows labeled by trt_factor (e.g., "Control" and "Treatment")
  effect_diff <- ale_data$.value[ale_data$trt_factor == "Treatment"] -
    ale_data$.value[ale_data$trt_factor == "Control"]
  
  return(effect_diff)
}

#######################################
# OUTER SIMULATION LOOP TO ASSESS COVERAGE
#######################################

# Set the number of simulation loops
n_sim <- 100  # You can change this to the desired number of simulation replications

# Vectors to store whether the 95% interval contains the true effect (0.15) in each simulation
coverage_weighted <- numeric(n_sim)
coverage_unweighted <- numeric(n_sim)

# Start simulation loop
for (sim in 1:n_sim) {
  cat("\n------------------------------\n")
  cat("Starting simulation loop", sim, "\n")
  
  # -----------------------------------------------------------
  # (a) Generate new replicate weights using a new seed per sim
  # -----------------------------------------------------------
  set.seed(1234 + sim)  # Set a new seed for each simulation iteration
  rep_design <- as.svrepdesign(design_obj, type = "bootstrap", replicates = 50)
  repweights <- weights(rep_design, type = "analysis")
  
  # Number of replicate weights (should be 50)
  B <- ncol(repweights)
  
  # Initialize vectors to store ALE differences for each replicate
  ale_diff_weighted <- numeric(B)
  ale_diff_unweighted <- numeric(B)
  
  # -----------------------------------------------------------
  # (b) Fit Weighted Random Forest Models with Replicate Weights
  # -----------------------------------------------------------
  for (i in 1:B) {
    cat("Fitting weighted model for replicate", i, "of simulation", sim, "\n")
    
    # Set seed (adding sim to ensure variability across simulations)
    set.seed(100 + i + sim)
    model_weighted <- train(formula_rf, 
                            data = survey_data, 
                            method = "ranger", 
                            weights = repweights[, i],
                            trControl = cv_control,
                            replace = TRUE, 
                            num.trees = 50,         
                            tuneGrid = tune_grid,
                            metric = "Specificity")
    
    # Compute the ALE difference for the weighted model
    ale_diff_weighted[i] <- compute_ALE_diff(model_weighted, survey_data)
  }
  
  # -----------------------------------------------------------
  # (c) Fit Unweighted Random Forest Models (Ignoring weights)
  # -----------------------------------------------------------
  for (i in 1:B) {
    cat("Fitting unweighted model for replicate", i, "of simulation", sim, "\n")
    
    set.seed(200 + i + sim)
    model_unweighted <- train(formula_rf, 
                              data = survey_data, 
                              method = "rf", 
                              importance = TRUE,
                              weights = rep(1, nrow(survey_data)),  # constant weights = unweighted
                              trControl = cv_control,
                              ntree = 50,
                              nodesize = 10,
                              tuneLength = 3)
    
    # Compute the ALE difference for the unweighted model
    ale_diff_unweighted[i] <- compute_ALE_diff(model_unweighted, survey_data)
  }
  
  # -----------------------------------------------------------
  # (d) Compute 95% Credible Intervals for Each Method
  # -----------------------------------------------------------
  lower_weighted <- quantile(ale_diff_weighted, 0.025)
  upper_weighted <- quantile(ale_diff_weighted, 0.975)
  
  lower_unweighted <- quantile(ale_diff_unweighted, 0.025)
  upper_unweighted <- quantile(ale_diff_unweighted, 0.975)
  
  # -----------------------------------------------------------
  # (e) Check if the True Effect (0.15) is Covered by Each Interval
  # -----------------------------------------------------------
  coverage_weighted[sim] <- ifelse(lower_weighted <= 0.15 && 0.15 <= upper_weighted, 1, 0)
  coverage_unweighted[sim] <- ifelse(lower_unweighted <= 0.15 && 0.15 <= upper_unweighted, 1, 0)
  
  cat("Simulation", sim, "Weighted interval: [", round(lower_weighted, 3), ",", round(upper_weighted, 3), "] -> Contains 0.15:", coverage_weighted[sim], "\n")
  cat("Simulation", sim, "Unweighted interval: [", round(lower_unweighted, 3), ",", round(upper_unweighted, 3), "] -> Contains 0.15:", coverage_unweighted[sim], "\n")
}

# -----------------------------------------------------------
# (f) Compute and Report Overall Coverage Probabilities
# -----------------------------------------------------------
coverage_prob_weighted <- mean(coverage_weighted)
coverage_prob_unweighted <- mean(coverage_unweighted)

cat("\n=====================================\n")
cat("Overall Coverage Probability (Weighted RF):", coverage_prob_weighted, "\n")
cat("Overall Coverage Probability (Unweighted RF):", coverage_prob_unweighted, "\n")
cat("=====================================\n")
