library(survival)
library(readxl)
library(ggRandomForests)
library(survminer)
library(zoo)
library(here)
library(caret)
library(DMwR2)
library(RANN)
library(e1071)
library(ranger)
library(tidyverse)
library(mice)
library(ROSE)
library(survey)

responses <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")
predictors <- c("DEPR_TOT", "AGE_BIN", "RIAGENDR", "RIDRETH1", "SMQ020",
                 "HIQ011", "BPQ020", "BMXBMI", "BMI_LVL", "BPQ080",
                "MEDDEP", "DMDBORNT", "PAQMV", "CADTOT", "DIDTOT",
                 "INC_BIN", "INC3", "DMDEDUC2", "CDQ008", "CDQ010", "TOT_REG")

# Combine responses and predictors
variables <- c(responses, predictors, "SDDSRVYR", "SDMVPSU", "SDMVSTRA", "MEC15YR", "SEQN")


# Subset the dataframe
subset_data <- subsetd_2 %>% select(all_of(variables))

#impute missing data using KNN method - ensure that columns are properly classified 
complete_data <- knnImputation(subset_data)
str(complete_data)
complete_data <- as.data.frame(complete_data)


# Survey design setup
design <- svydesign(
  id = ~SDMVPSU,
  weights = ~MEC15YR,
  strata = ~SDMVSTRA,
  data = complete_data,
  nest = TRUE
)

set.seed(123)
rep_design <- as.svrepdesign(design, type = "bootstrap", replicates = 10)
rep_design <- subset(rep_design, !apply(
  rep_design$variables[, responses], 1,
  function(x) all(x == 0 | is.na(x))
))

# Stabilize weights
repweights <- weights(rep_design, type = "analysis") / mean(weights(rep_design, type = "analysis"))
pweights <- weights(rep_design, type = "sampling") / mean(weights(rep_design, type = "sampling"))

pweights <- as.numeric(pweights)
data.y <- rep_design$variables[,c(responses)]
data <- rep_design$variables[,c(responses, predictors)]


survival_forest <- rfsrc(get.mv.formula(responses), 
                         data = data, 
                         ntree = 100,
                         save.memory = FALSE,
                         nsplit = 10,
                         importance = TRUE,
                         case.wt = pweights,
                         splitrule = "mahalanobis")