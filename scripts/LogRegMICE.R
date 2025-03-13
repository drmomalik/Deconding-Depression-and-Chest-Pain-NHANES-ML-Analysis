library(survey)
library(mice)
library(howManyImputations)

responses <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")
predictors <- c("DEPR_BIN", "DEPR_LVL", "AGE_BIN", "RIAGENDR", "RIDRETH1", "SMQ020",
                "SMQ040", "HIQ011", "BPQ020", "BMXBMI", "BMI_LVL", "BPQ080",
                "ALQ130", "MEDDEP", "DMDBORNT", "PAQMV", "CADTOT", "DIDTOT",
                "DUQTOT", "INC_BIN")

# subset to only include the variables we are working with to improve computational efficiency
wd_subset <- wd[, c(responses, predictors, "SDMVPSU", "SDMVSTRA", "MEC15YR", "SEQN")]

# Initialize the predictor matrix
predictorMatrix <- mice::make.predictorMatrix(wd_subset)

# Set all columns except the predictors to 0 (do not use them to impute)
predictorMatrix[, !colnames(predictorMatrix) %in% predictors] <- 0

# Ensure only predictor variables are imputed
predictorMatrix[!rownames(predictorMatrix) %in% predictors, ] <- 0

# Run mice with the custom predictor matrix
imputations <- mice(wd_subset, seed = 123, maxit = 5, m = 5, predictorMatrix = predictorMatrix,
                    printFlag = TRUE)


lm_svy_mi <- function(formula, imputations) {
  
  #setting up null objects allows us to easily add results
  #later
  b <- se <- NULL
  
  #now loop through our imputations and run the model
  for(i in 1:imputations$m) {
    #grab the complete dataset
    imputation <- complete(imputations, i)
    # Normalize weights so they sum to the total sample size
    imputation$MEC15YR <- imputation$MEC15YR / sum(imputation$MEC15YR, na.rm = TRUE) * nrow(imputation)
    # Verify normalization
    sum(imputation$MEC15YR, na.rm = TRUE)
    #create the design effect object
    imputation.svy <- svydesign(
      id = ~SDMVPSU, 
      weights = ~MEC15YR, 
      strata = ~SDMVSTRA, 
      data = imputation,
      nest = TRUE
    ) 
    
    subset_data <- subset(
      imputation.svy,  # Access the data (variables) in the svydesign object
      !apply(
        imputation.svy$variables[, c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")],
        1, 
        function(x) all(x == 0 | is.na(x))  # Keep rows where not all values of CDQ009A-CDQ009H are 0 or NA
      )
    )
    
    #run the model
    model <- svyglm(formula, design=subset_data, family = quasibinomial)
    #collect the results
    b <- cbind(b, coef(model))
    se <- cbind(se, summary(model)$coef[,2])

  }
  
  #now pool the results
  b.pool <- apply(b, 1, mean)
  or <- exp(b.pool)
  between.var <- apply(b, 1, function(x) var(x))
  within.var <- apply(se^2, 1, mean)
  se.pool <- sqrt(within.var+between.var+between.var/imputations$m) 
  t.pool <- b.pool/se.pool 
  pvalue.pool <- (1-pnorm(abs(t.pool)))*2 
  coefficients <- data.frame(b.pool, se.pool, t.pool, pvalue.pool)
  lci <- exp(b.pool-1.96*se.pool)
  uci <- exp(b.pool+1.96*se.pool)
  
  
  #return everything in a list
  return(list(coef=coefficients,
              OR=or,
              Lower.CI=lci,
              Upper.CI=uci,
              n=n))
}



#Loop over all responses 
models <- list()

for (response in responses) {
  formula <- as.formula(paste(response, " ~ ", paste(predictors, collapse = " + ")))
  models[response] <- lm_svy_mi(formula, imputations)
}
