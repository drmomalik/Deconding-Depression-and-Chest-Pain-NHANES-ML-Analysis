
library(glmnet)

# Fit separate penalized logistic regression models for each outcome
models <- lapply(
  c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", 
    "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H"),
  function(outcome) {
    x <- model.matrix(~ DEPR_BIN + CADTOT + MEDDEP, data)[, -1]
    y <- data[[outcome]]
    glmnet(x, y, alpha = 0, lambda = 0, family = "binomial")
  }
)

predicted_probs <- sapply(models, function(model) {
  x <- model.matrix(~ DEPR_BIN + CADTOT + MEDDEP, data)[, -1]
  predict(model, newx = x, type = "response")
})
colnames(predicted_probs) <- c("Right_Chest", "Right_Arm", "Neck", "Left_Chest",
                               "Left_Arm", "Jaw", "Back", "Epigastrium")

# Correlation matrix of predicted probabilities
cor_pred_probs <- cor(predicted_probs)

# Alternatively, correlation of observed binary outcomes
observed_outcomes <- data[, c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D",
                              "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")]
cor_observed <- cor(observed_outcomes)

library(corrplot)

# Visualize correlation of predicted probabilities
corrplot(cor_pred_probs, method = "color", type = "upper", tl.cex = 0.8,
         title = "Correlation of Predicted Probabilities")

# Visualize correlation of observed outcomes
corrplot(cor_observed, method = "color", type = "upper", tl.cex = 0.8,
         title = "Correlation of Observed Outcomes")


# Subset data by depression status
data_depressed <- data[data$DEPR_BIN == 1, ]
data_non_depressed <- data[data$DEPR_BIN == 0, ]

# Correlation matrices for each group
cor_depressed <- cor(data_depressed[, c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D",
                                        "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")])
cor_non_depressed <- cor(data_non_depressed[, c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D",
                                                "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")])

# Visualize side by side
par(mfrow = c(1, 2))
corrplot(cor_depressed, method = "color", type = "upper", tl.cex = 0.8, title = "Depressed")
corrplot(cor_non_depressed, method = "color", type = "upper", tl.cex = 0.8, title = "Non-Depressed")


library(vegan)

# Mantel's test for correlation matrices
mantel_test <- mantel(cor_depressed, cor_non_depressed, method = "pearson", permutations = 999)
print(mantel_test)

