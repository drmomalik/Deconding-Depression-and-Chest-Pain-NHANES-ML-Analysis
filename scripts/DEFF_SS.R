## Sample size estimate

## Assume a conservative deff = 2
library(pwr)
pwr.t.test(n = NULL, d = 0.3, sig.level = 0.05, power = 0.8, 
           type = c("two.sample"),
           alternative = c("two.sided"))

## Calculate DEFF for CAD and other variables of interest
# Define variables of interest
variables <- c("CAD", "DEPR_NEW", "CDQ009A", "CDQ009B", "CDQ009C", 
               "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")

# Initialize an empty data frame to store results
deff_results <- data.frame(Variable = character(), DEFF = numeric(), stringsAsFactors = FALSE)

# Loop through each variable and compute DEFF
for (var in variables) {
  # Compute survey mean with DEFF
  svy_mean <- svymean(as.formula(paste("~", var)), matched_design, deff = TRUE)
  
  # Extract DEFF and store results
  deff_value <- attr(svy_mean, "deff")[1]  # Extract the DEFF for the variable
  deff_results <- rbind(deff_results, data.frame(Variable = var, DEFF = deff_value))
}

# Print results
print(deff_results)
  


  