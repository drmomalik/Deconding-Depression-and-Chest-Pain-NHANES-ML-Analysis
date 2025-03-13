library(randomForest)
library(randomForestSRC)
library(mice)
library(survey)
library(pROC)


responses <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")
predictors <- c("DEPR_TOT","AGE_BIN", "RIAGENDR", "RIDRETH1", "SMQ020",
                "SMQ040", "HIQ011", "BPQ020", "BMXBMI", "BMI_LVL", "BPQ080",
                "ALQ130", "MEDDEP", "DMDBORNT", "PAQMV", "CADTOT", "DIDTOT",
                "DUQTOT", "INC_BIN", "INC3", "DMDEDUC2", "CDQ008", "CDQ010")
wd_subset <- wd[, c(responses, predictors,"SDDSRVYR", "SDMVPSU", "SDMVSTRA", "MEC15YR", "SEQN")]


# Survey design setup
design <- svydesign(
  id = ~SDMVPSU,
  weights = ~MEC15YR,
  strata = ~SDMVSTRA,
  data = wd_subset,
  nest = TRUE
)

set.seed(123)
rep_design <- as.svrepdesign(design, type = "bootstrap", replicates = 500)
rep_design <- subset(rep_design, !apply(
  rep_design$variables[, responses], 1,
  function(x) all(x == 0 | is.na(x))
))

# Stabilize weights
repweights <- weights(rep_design, type = "analysis") / mean(weights(rep_design, type = "analysis"))
pweights <- weights(rep_design, type = "sampling") / mean(weights(rep_design, type = "sampling"))


data <- rep_design$variables[c(predictors, responses)]


# Define the function for processing each replicate
process_replicate <- function(i, repweights, data, responses, predictors) {
  
  vimp_rep <- matrix(0, ncol = length(responses), nrow = length(predictors)+7)
  perf_rep <- matrix(0, ncol = length(responses), nrow = 20)
  
  current_weights <- repweights[,i]
  
  for (j in 1:length(responses)) {
    formula <- as.formula(paste(responses[j], "~."))
    
    
    # Use of imbalanced quantile regressor, higher tree used, and case weights applied to account for survey design
    rf_model <- imbalanced.rfsrc(formula, data = data, importance = "permute",
                                 samptype = "swr", case.wt = current_weights, na.action = "na.impute", 
                                 nimpute = 5, n = 2000)
    
    vimp_rep[, j] <- rf_model$importance[, 1]
    perf_rep[, j] <- get.imbalanced.performance(rf_model)
  }
  
  return(list(vimp_rep = vimp_rep, perf_rep = perf_rep))
}

# Initialize storage for results
vimp_results <- list()
perf_results <- list()

# Loop through each replicate
for (i in 1:ncol(repweights)) {
  cat(sprintf("Processing replicate %d of %d...\n", i, ncol(repweights)))  # Track progress
  
  # Call the processing function
  replicate_result <- process_replicate(i, repweights, data, responses, predictors)
  
  # Store results
  vimp_results[[i]] <- replicate_result$vimp_rep
  perf_results[[i]] <- replicate_result$perf_rep
}

cat("All replicates processed.\n")


# Aggregate the results
vimp_list <- lapply(vimp_results, function(res) res$vimp_rep)
perf_list <- lapply(perf_results, function(res) res$perf_rep)

# Convert to matrices
vimp_matrix <- do.call(cbind, vimp_list)
perf_matrix <- do.call(rbind, perf_list)

# Aggregate (mean, sd)
vimp_mean <- rowMeans(vimp_matrix)
vimp_sd <- apply(vimp_matrix, 1, sd)

perf_mean <- rowMeans(perf_matrix)
perf_sd <- apply(perf_matrix, 1, sd)

# Print results
print(vimp_mean)
print(vimp_sd)
print(perf_mean)
print(perf_sd)


