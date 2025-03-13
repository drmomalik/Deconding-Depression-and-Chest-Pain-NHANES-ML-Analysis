library(rpms)
library(MLPUGS)
library(survey)
library(mice)
library(randomForest)
library(iml)
library(tidyverse)


responses <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")
predictors <- c("DEPR_TOT","AGE_BIN", "RIAGENDR", "RIDRETH1", "SMQ020",
                "HIQ011", "BPQ020", "BMXBMI", "BPQ080",
                "MEDDEP", "DMDBORNT", "PAQMV", "CADTOT", "DIDTOT",
                "INC3", "DMDEDUC2", "CDQ008", "CDQ010")
wd_subset <- wd[, c(responses, predictors,"SDDSRVYR", "SDMVPSU", "SDMVSTRA", "MEC15YR", "SEQN")]


# Initialize the predictor matrix for MICE
predictorMatrix <- mice::make.predictorMatrix(wd_subset)
predictorMatrix[, !colnames(predictorMatrix) %in% predictors] <- 
predictorMatrix[!rownames(predictorMatrix) %in% predictors, ] <- 0

# Run MICE
imputations <- mice(wd_subset, seed = 123, maxit = 1, m = 1, predictorMatrix = predictorMatrix, printFlag = TRUE)

# Run RF imputations (missing Forest)



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



data <- complete(imputations, 1)
complete_data <- complete(imputations, 1)
data <- subset(complete_data, !apply(
  complete_data[, responses], 1,
  function(x) all(x == 0 | is.na(x))
))
data <- data[,c(predictors, responses)]
df <- data
x <- data[, c(predictors)]
y <- data[, c(responses)]
weights <- data$MEC15YR
strata <- data$SDMVSTRA
newdata <- x


model <- ecc(x, y, m=3, silent=FALSE, .f = randomForest::randomForest, weights=weights, strata=strata)
predictions <- predict(model, newdata = x, n.iters = 1000, burn.in = 100, silent=FALSE,
              .f = function(rF, newdata){ randomForest:::predict.randomForest(rF, type="prob")})
summary(predictions)

y_num <- as.data.frame(lapply(y, function(x) {as.numeric(as.character(x))}))

baseline_perf <- validate_pugs(predictions, y_num)
hamming_baseline <- as.numeric(baseline_perf)[5]

hamming_loss <- function(predictions, y_num) {
  all_scores <- validate_pugs(predictions, y_num)
  hamming_loss <- as.numeric(all_scores)[5]
  
  return(hamming_loss)
}





## VIMP
# Function to compute joint VIMP by permuting each predictor
compute_joint_permutation_importance <- function(model, data, y_num, predictors, hamming_baseline, n_permutations = 10) {
  importance_scores <- numeric(length(predictors))  # To store the aggregated importance for each predictor
  
  for (i in 1:length(predictors)) {
    predictor <- predictors[i]
    
    # Initialize a vector to store performance drops for multiple permutations
    performance_drops <- numeric(n_permutations)
    
    # Loop through multiple permutations for each predictor
    for (j in 1:n_permutations) {
      print(paste("Starting permutation #",j,"for",predictor))
      # Save the original values of the predictor column
      original_values <- data[[predictor]]
      
      # Permute the predictor (shuffle its values)
      data[[predictor]] <- sample(data[[predictor]])
      
      # Calculate performance after permutation (for the joint model)
      predictions_permuted <- predict(model, newdata = data, n.iters = 100, burn.in = 10,thin=2,silent=FALSE,
                                      .f = function(rF, newdata){ randomForest:::predict.randomForest(rF, type="prob")})
      
      # Define Hamming loss function
      hamming_loss <- function(predictions, y_num) {
        all_scores <- validate_pugs(predictions, y_num)
        hamming_loss <- as.numeric(all_scores)[5]
        return(hamming_loss)
      }
      
      # Calculate the permuted performance
      permuted_performance <- hamming_loss(predictions_permuted, y_num)
      
      # Calculate the performance drop for this permutation
      performance_drop <- permuted_performance - hamming_baseline
      performance_drops[j] <- performance_drop
      
      # Restore the original values of the predictor
      data[[predictor]] <- original_values
    }
    
    # Aggregate the performance drops for the predictor (e.g., by calculating the mean)
    importance_scores[i] <- mean(performance_drops)
    names(importance_scores)[i] <- predictor
    print(paste("VIMP for",predictor,":",importance_scores[i]))
  }
  
  # Return the aggregated importance scores for each predictor
  return(importance_scores)
}


# List of predictors to evaluate (adjust based on your dataset)
predictors <- c("DEPR_TOT","AGE_BIN", "RIAGENDR", "RIDRETH1", "SMQ020",
                "HIQ011", "BPQ020", "BMXBMI", "BPQ080",
                "MEDDEP", "DMDBORNT", "PAQMV", "CADTOT", "DIDTOT",
                "INC3", "DMDEDUC2", "CDQ008", "CDQ010")

# Calculate joint permutation importance
vimp_joint <- compute_joint_permutation_importance(model=model, data=x, y_num=y_num, 
            predictors=predictors, hamming_baseline=hamming_baseline,
            n_permutations=50)

# Store results
vimp_joint


