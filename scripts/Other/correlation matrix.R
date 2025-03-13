
dfc <- wd

## Correlation matrix for regions of chest pain stratified by PHQ9-score

# List of binary variables
binary_vars <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")

# Convert binary variables to numeric
dfc[binary_vars] <- lapply(dfc[binary_vars], as.numeric)

# Function to calculate correlation matrix for a subset
calculate_correlation <- function(subset_df) {
  cor(subset_df[, binary_vars], method = "pearson")  # Can use "spearman" if desired
}

# Split the data by the categorical variable
stratified_correlation <- lapply(split(dfc, dfc$DEPR_LVL), calculate_correlation)

# View correlation matrices for each group
stratified_correlation

library(corrplot)

# Desired labels for the correlation matrix
labels <- c(
  "Right Arm", 
  "Right Chest", 
  "Neck", 
  "Upper Sternum", 
  "Lower Sternum", 
  "Left Chest", 
  "Left Arm", 
  "Epigastric"
)

# Visualize correlation matrix for each level with custom labels
for (level in names(stratified_correlation)) {
  cat("Correlation Matrix for Group:", level, "\n")
  
  corr_matrix <- stratified_correlation[[level]]
  
  # Set the custom row and column labels
  colnames(corr_matrix) <- labels
  rownames(corr_matrix) <- labels
  
  # Plot the correlation matrix
  corrplot(
    corr_matrix, 
    method = "circle", 
    tl.col = "black",  # Text color
    tl.cex = 0.8,      # Text size
    mar = c(0, 0, 2, 0) # Adjust margin for the title
  )
  
  # Add the title manually
  title(
    main = paste("Depr Lvl", level),
    line = -2,  # Moves the title closer to the plot (negative values bring it closer)
    cex.main = 1.2 # Adjust title font size if needed
  )
}



## Correlation matrix for PHQ-9 score (cont.) by total regions of chest pain

# Calculate correlation while excluding NAs
correlation <- cor(wd$TOT_REG, wd$DEPR_TOT, use = "complete.obs", method = "pearson")
print(correlation)

library(ggplot2)

# Create a data frame excluding rows with NAs
wd_clean <- na.omit(wd[, c("TOT_REG", "DEPR_TOT")])

# Scatterplot with regression line
ggplot(wd_clean, aes(x = TOT_REG, y = DEPR_TOT)) +
  geom_point(color = "blue", size = 2) +
  geom_smooth(method = "lm", color = "red") +  # Add regression line
  labs(title = paste("Correlation =", round(correlation, 2)),
       x = "TOT_REG", y = "DEPR_TOT") +
  theme_minimal()



### survey weighted corr matrix 
library(survey)
library(jtools)
library(corrplot)

dfc <- wd

# List of binary variables
binary_vars <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")

# Convert binary variables to numeric
dfc[, binary_vars] <- lapply(dfc[, binary_vars], as.numeric)

# Define the survey design
sub_design <- svydesign(id = ~SDMVPSU, weights = ~MEC15YR,
                        strata = ~SDMVSTRA, data = dfc, nest = TRUE)

# Subset the survey design based on conditions
sub_design <- subset(sub_design, !apply(
  sub_design$variables[, c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")],
  1, 
  function(x) all(x == 0 | is.na(x))  # Keep rows where not all values of CDQ009A-CDQ009H are 0 or NA
))

# Create two subsets based on DEPR_TOT
des_1 <- subset(sub_design, DEPR_TOT < 15)
des_2 <- subset(sub_design, DEPR_TOT >= 15)

# List of designs to loop over
designs <- list("PHQ-9 Score < 15" = des_1, "PHQ-9 Score ≥ 15" = des_2)

# Desired labels for the correlation matrix
labels <- c(
  "Right Arm", 
  "Right Chest", 
  "Neck", 
  "Upper Sternum", 
  "Lower Sternum", 
  "Left Chest", 
  "Left Arm", 
  "Epigastric"
)

# Loop through both designs
for (level in names(designs)) {
  # Get the design for the current level
  current_design <- designs[[level]]
  
  # Calculate survey-weighted correlation for the current design
  corr_result <- svycor(
    ~CDQ009A + CDQ009B + CDQ009C + CDQ009D + CDQ009E + CDQ009F + CDQ009G + CDQ009H,
    design = current_design,
    na.rm = TRUE,
    sig.stats = TRUE
  )
  
  # Extract the correlation matrix from the svycor object
  corr_matrix <- corr_result$cor
  
  # Set custom row and column labels
  colnames(corr_matrix) <- labels
  rownames(corr_matrix) <- labels
  
  # Plot the correlation matrix for the current group
  corrplot(
    corr_matrix, 
    method = "circle", 
    tl.col = "black",  # Text color
    tl.cex = 0.8,      # Text size
    mar = c(0, 0, 2, 0) # Adjust margin for the title
  )
  
  # Add the title manually
  title(
    main = paste("Correlation Matrix -", level),  # Add group-specific title
    line = 2.5,  # Moves the title closer to the plot (negative values bring it closer)
    cex.main = 1.1 # Adjust title font size if needed
  )
}


### correlation b/w DEPR_TOT and TOT_REG

library(survey)
library(jtools)
library(corrplot)

dfc <- wd


# Define the survey design
sub_design <- svydesign(id = ~SDMVPSU, weights = ~MEC15YR,
                        strata = ~SDMVSTRA, data = dfc, nest = TRUE)

# Calculate survey-weighted correlation
corr_result <- svycor(
  ~TOT_REG + DEPR_TOT,
  design = sub_design,
  na.rm = TRUE,
  sig.stats = TRUE,
  digits = 6
)

# Extract the correlation matrix from the svycor object
corr_matrix <- corr_result$cor

# Desired labels for the correlation matrix
labels <- c(
  "# of regions of pain",
  "PHQ-9 SCORE"
)

# Set custom row and column labels
colnames(corr_matrix) <- labels
rownames(corr_matrix) <- labels

# Plot the correlation matrix
corrplot(
  corr_matrix, 
  method = "circle", 
  tl.col = "black",  # Text color
  tl.cex = 0.8,      # Text size
  mar = c(0, 0, 2, 0) # Adjust margin for the title
)

# Add the title manually
title(
  main = "Correlation Matrix - PHQ-9 score and # of regions of pain",  # Adjust title as needed
  line = -2,  # Moves the title closer to the plot (negative values bring it closer)
  cex.main = 1.2 # Adjust title font size if needed
)

