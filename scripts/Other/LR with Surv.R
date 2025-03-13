# Load necessary libraries
library(survey)

# Create the survey design object
svydesign_object <- svydesign(
  id = ~SDMVPSU, 
  weights = ~MEC15YR, 
  strata = ~SDMVSTRA, 
  data = wd,
  nest = TRUE
)

# Subset the design object based on the conditions
svydesign_object <- subset(
  svydesign_object, 
  !is.na(CDQ001) & CDQ001 == 1
)

subset_data <- subset(
  svydesign_object,  # Access the data (variables) in the svydesign object
  !apply(
    svydesign_object$variables[, c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")],
    1, 
    function(x) all(x == 0 | is.na(x))  # Keep rows where not all values of CDQ009A-CDQ009H are 0 or NA
  )
)


# Define the variables for the logistic regression
responses <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")
predictors <- c("DEPR_LVL", "AGE_BIN", "RIAGENDR", "RIDRETH1", 
                "SMQ020", "SMQ040", "HIQ011", "BPQ020", "BMXBMI", "BPQ080",
                "ALQ130", "MEDDEP", "DMDBORNT", "PAQMV", "CADTOT", "DIDTOT",
                "DUQTOT", "INC_BIN")

# Initialize an empty list to store the regression models
regression_models <- list()

# Run logistic regression for each response
for (response in responses) {
  formula <- as.formula(paste(response, "~ ", paste(predictors, collapse = " + ")))
  regression_models[[response]] <- svyglm(formula, design = subset_data, family = quasibinomial())
}

# Extract residuals from each model
residuals_list <- lapply(regression_models, residuals)

# Create a correlation matrix of the residuals
residuals_matrix <- do.call(cbind, residuals_list)
correlation_matrix <- cor(residuals_matrix, use = "complete.obs")

# Print the correlation matrix
print(correlation_matrix)