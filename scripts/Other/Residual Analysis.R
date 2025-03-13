#Residuals vs. Fitted Values:

#  Purpose: Checks for patterns indicating non-linearity or heteroscedasticity.
#Expectation: The residuals should be randomly scattered around zero.

#Residuals vs. Predictors:

##  Purpose: Examines if a specific predictor is inadequately modeled.
#Expectation: No discernible pattern in the residuals.

#Normal Q-Q Plot of Residuals:

#  Purpose: Checks the normality of residuals.
#Expectation: Residuals should align with the diagonal line for a well-fitting model.

#Scale-Location Plot (Spread-Location Plot):

#  Purpose: Tests for homoscedasticity (constant variance).
#Expectation: The spread of residuals should remain roughly constant across fitted values.

#Leverage vs. Residuals (Cook's Distance):

#Purpose: Identifies influential observations.
#Expectation: No single observation should dominate the model.


# Set variance approximation for single-PSU strata
options(survey.lonely.psu = "adjust")

# Set up the survey design
nhanes_design <- svydesign(
  id = ~SDMVPSU,         # Primary sampling units
  strata = ~SDMVSTRA,    # Sampling strata
  weights = ~MEC15YR,   # MEC weights
  data = wd,    # Your NHANES dataframe
  nest = TRUE            # Indicates nesting of PSU within strata
)

# Ex: Logistic regression: Predict CDQ009A using DEPR_LVL
model0 <- svyglm(CDQ009A ~ DEPR_LVL, 
                 design = nhanes_design, 
                 family = quasibinomial())  # Use quasibinomial for survey-adjusted logit
model1 <- svyglm(CDQ009A ~ DEPR_LVL+BMXBMI+INDFMPIR+RIDAGEYR+ALQ130+
                  DIDTOT+CAD+PAQMV+DMDBORNT+MEDDEP+SDDSRVYR+BPQ080+
                  BPQ020+HIQ011+SMQ020+RIDRETH1+RIAGENDR, 
                design = nhanes_design, 
                family = quasibinomial())  # Use quasibinomial for survey-adjusted logit
anova(model1, test="F")
summary(model1)
residuals(model1)


model <- model1

# Residuals vs. Fitted
plot(model$fitted.values, residuals(model, type = "deviance"),
     main = "Residuals vs Fitted",
     xlab = "Fitted values",
     ylab = "Deviance Residuals")
abline(h = 0, col = "red", lwd = 2)

# Normal Q-Q Plot
qqnorm(residuals(model, type = "deviance"),
       main = "Normal Q-Q Plot")
qqline(residuals(model, type = "deviance"), col = "red")

# Scale-Location Plot
sqrt_abs_residuals <- sqrt(abs(residuals(model, type = "deviance")))
plot(model$fitted.values, sqrt_abs_residuals,
     main = "Scale-Location",
     xlab = "Fitted values",
     ylab = "√|Residuals|")
abline(h = 0, col = "red")



### Residuals for each predictor

dataset <- wd

# Extract deviance residuals
residuals_deviance <- residuals(model, type = "deviance")

# List of predictors
predictors <- names(dataset)[names(dataset) %in% c("DEPR_LVL","BMXBMI","INDFMPIR","RIDAGEYR","ALQ130",
                  "DIDTOT","CAD","PAQMV","DMDBORNT","MEDDEP","SDDSRVYR","BPQ080",
                  "BPQ020","HIQ011","SMQ020","RIDRETH1","RIAGENDR")]

# Loop through predictors
for (predictor in predictors) {
  # Check if predictor is numeric or categorical
  if (is.numeric(dataset[[predictor]])) {
    # Scatter plot for numeric predictors
    plot(dataset[[predictor]], residuals_deviance,
         main = paste("Residuals vs", predictor),
         xlab = predictor,
         ylab = "Deviance Residuals",
         pch = 19, col = "blue")
    abline(h = 0, col = "red", lwd = 2)
  } else {
    # Boxplot for categorical predictors
    boxplot(residuals_deviance ~ dataset[[predictor]],
            main = paste("Residuals vs", predictor),
            xlab = predictor,
            ylab = "Deviance Residuals",
            col = "lightblue")
    abline(h = 0, col = "red", lwd = 2)
  }
}




