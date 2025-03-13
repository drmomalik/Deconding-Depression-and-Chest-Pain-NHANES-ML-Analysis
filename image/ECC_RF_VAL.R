library(survey)
library(splitTools)

# Define outcome variables
outcome_cols <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", 
                  "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")

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


train_proportion <- 0.67
set.seed(123)

# Split data by strata
split_data <- lapply(split(survey_design$variables, survey_design$variables$SDMVSTRA), function(stratum_data) {
  train_indices <- sample(seq_len(nrow(stratum_data)), size = floor(train_proportion * nrow(stratum_data)), replace = FALSE)
  list(train = stratum_data[train_indices, ], test = stratum_data[-train_indices, ])
})

# Combine the train and test datasets from each stratum
train_data <- do.call(rbind, lapply(split_data, function(x) x$train))
test_data <- do.call(rbind, lapply(split_data, function(x) x$test))

# Create the survey design objects
train_design <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~MEC15YR, data = train_data, nest=TRUE)
test_design <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~MEC15YR, data = test_data, nest=TRUE)

# Create a train rep and test rep designs for validation of ECC model 
train_repwt <- as.svrepdesign(train_design, type="bootstrap", replicates=100)
test_rep <- as.svrepdesign(test_design, type="bootstrap", replicates=100)

# Need to create actual weights in order to bootstrap while training (do not need to do for testing)
train_repweights <- weights(train_repwt, type = "analysis") / mean(weights(train_repwt, type = "analysis"))


