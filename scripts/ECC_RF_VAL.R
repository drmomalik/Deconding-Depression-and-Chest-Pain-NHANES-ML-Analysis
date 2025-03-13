library(survey)
library(missForest)
library(splitTools)
library(rsample)

# Define outcome variables
responses <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")

predictors <- c("DEPR_TOT","AGE_BIN", "RIAGENDR", "RIDRETH1", "SMQ020",
                "SMQ040", "HIQ011", "BPQ020", "BMXBMI", "BPQ080",
                "ALQ130", "MEDDEP", "DMDBORNT", "PAQMV", "CADTOT", "DIDTOT",
                "DUQTOT", "INC3", "DMDEDUC2", "CDQ008", "CDQ010")
wd_subset <- wd[, c(responses, predictors,"SDDSRVYR", "SDMVPSU", "SDMVSTRA", "MEC15YR", "SEQN")]


# Extract the data from the svrepdesign object
design <- svydesign(
  id = ~SDMVPSU,
  weights = ~MEC15YR,
  strata = ~SDMVSTRA,
  data = wd_subset,
  nest = TRUE
)

survey_design <- subset(design, !apply(
  design$variables[, responses], 1,
  function(x) all(x == 0 | is.na(x))
))

survey_design$variables$SDMVSTRA <- NULL
survey_design$variables$SDMVPSU <- NULL
survey_design$variables$SEQN <- NULL

imputed_data <- missForest(xmis=survey_design$variables, verbose = TRUE)
# Extract the completed dataset
wd_subset_imputed <- imputed_data$ximp  # This is the dataset with imputed values
print(imputed_data$OOBerror)

survey_design$variables <- data <- wd_subset_imputed

split <- initial_split(survey_design$variables, strata = "CDQ009H", prop = 0.75)
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

current_label <- "CADTOT"
formula <- as.formula(paste(current_label, "~", paste(c(predictors, responses), collapse="+")))
data <- train_data

data$weights <- train_repweights[,1]

data_minority <- data %>% filter(CADTOT == 1)
data_majority <- data %>% filter(CADTOT == 0)
ratio <- nrow(data_majority)/nrow(data_minority)

beta <- 2
majority_sampled <- data_majority %>% sample_n(400, replace = TRUE, weight=weights)
minority_sampled <- data_minority %>% sample_n(400, replace = TRUE, weight=weights)
balanced_data <- bind_rows(data_minority, majority_sampled)

model <- imbalanced(formula, data=balanced_data, case.weights = data$weights, 
                                      method="rfq", splitrule = "auc", ntree =3000, 
                                      nodesize=1, fast=FALSE, forest=TRUE, importance=TRUE)
pred<-predict(model, newdata=test_data, outcome="test")
th <- get.imbalanced.optimize(model)["threshold"]
get.imbalanced.performance(pred, threshold = th)
print(model)





model <- rfsrc(formula, data=balanced_data, case.weights = data$weights, 
                splitrule = "gini", ntree = 3000, 
               nodesize=1, fast=FALSE, forest=TRUE)
print(model)
model <- ranger(formula, data=balanced_data, num.trees = 1000, splitrule="gini",
        replace=TRUE, min.node.size = 1, probability = FALSE)
print(model)
print(model$confusion.matrix)

model <- glmnet(formula, test_design, family=quasibinomial())
summary(model)


library(randomForest)

current_label <- "CDQ009A"
data <- train_data
formula <- as.formula(paste(current_label, "~", paste(predictors, collapse="+")))

data$weights <- train_repweights[,1]

# Split data into classes
data_minority <- data %>% filter(CDQ009A == 1)
data_majority <- data %>% filter(CDQ009A == 0)
ratio <- nrow(data_majority)/nrow(data_minority)

beta <- sqrt(ratio)

# Initialize an empty randomForest object to store all trees
all <- NULL

# Loop to train 20 models and combine them
for (i in 1:50) {
  # Under-sample majority class to match minority size
  majority_sampled <- data_majority %>% sample_n(beta*nrow(data_minority), replace = TRUE)
  balanced_data <- bind_rows(data_minority, majority_sampled)
  
  # Train a randomForest model with case weights
  model <- randomForest(
    formula, 
    data = balanced_data, 
    ntree = 10,
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

true_label <- train_data$CDQ009A
predictions <- predict(all, newdata=train_data, type="response")

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