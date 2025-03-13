
## Difference distributions

# for depression cohort
# Get unique predictors
unique_predictors <- unique(ale_dep_sens$predictor)

# Initialize an empty list to store vectors
diff_list <- list()

# Loop through each predictor
for(pred in unique_predictors) {
  # Subset data for current predictor
  pred_data <- ale_dep_sens[ale_dep_sens$predictor == pred, ]
  
  # Get unique replicates
  replicates <- unique(pred_data$replicate)
  
  # Initialize vector for differences
  diffs <- numeric(length(replicates))
  
  # Calculate differences for each replicate
  for(i in seq_along(replicates)) {
    rep_data <- pred_data[pred_data$replicate == replicates[i], ]
    
    # Get values for level 1 and 2
    value1 <- rep_data$.value[rep_data$value == 1]
    value2 <- rep_data$.value[rep_data$value == 2]
    
    # Ensure values exist before subtraction
    if (length(value1) > 0 & length(value2) > 0) {
      diffs[i] <- value2 - value1
    } else {
      diffs[i] <- NA  # Handle missing values gracefully
    }
  }
  
  # Store the vector in a list
  diff_list[[pred]] <- diffs
}

# Convert list to data frame (ensuring columns align correctly)
diff_distributions_dep <- as.data.frame(diff_list)

#for no depression cohort

# Get unique predictors
unique_predictors <- unique(ale_nodep_sens$predictor)

# Initialize an empty list to store vectors
diff_list <- list()

# Loop through each predictor
for(pred in unique_predictors) {
  # Subset data for current predictor
  pred_data <- ale_nodep_sens[ale_nodep_sens$predictor == pred, ]
  
  # Get unique replicates
  replicates <- unique(pred_data$replicate)
  
  # Initialize vector for differences
  diffs <- numeric(length(replicates))
  
  # Calculate differences for each replicate
  for(i in seq_along(replicates)) {
    rep_data <- pred_data[pred_data$replicate == replicates[i], ]
    
    # Get values for level 1 and 2
    value1 <- rep_data$.value[rep_data$value == 1]
    value2 <- rep_data$.value[rep_data$value == 2]
    
    # Ensure values exist before subtraction
    if (length(value1) > 0 & length(value2) > 0) {
      diffs[i] <- value2 - value1
    } else {
      diffs[i] <- NA  # Handle missing values gracefully
    }
  }
  
  # Store the vector in a list
  diff_list[[pred]] <- diffs
}

# Convert list to data frame (ensuring columns align correctly)
diff_distributions_nodep <- as.data.frame(diff_list)

# Net effect 
library(ggplot2)
library(perm)

# Initialize storage
predictors <- colnames(diff_distributions_dep)
results <- data.frame(Predictor = predictors, 
                      p_value = numeric(length(predictors)), 
                      mean_diff = numeric(length(predictors)), 
                      p_value_adj = numeric(length(predictors)), 
                      ci_lower_dep = numeric(length(predictors)), 
                      ci_upper_dep = numeric(length(predictors)), 
                      mean_dep = numeric(length(predictors)), 
                      q1_dep = numeric(length(predictors)), 
                      q3_dep = numeric(length(predictors)), 
                      ci_lower_nodep = numeric(length(predictors)), 
                      ci_upper_nodep = numeric(length(predictors)), 
                      mean_nodep = numeric(length(predictors)), 
                      q1_nodep = numeric(length(predictors)), 
                      q3_nodep = numeric(length(predictors)))

# Loop over each predictor
for (i in seq_along(predictors)) {
  predictor <- predictors[i]
  print(paste("Processing predictor:", predictor))
  
  # Get bootstrap estimates
  dep_values <- diff_distributions_dep[[predictor]]
  nodep_values <- diff_distributions_nodep[[predictor]]
  
  # Permutation test using perm package
  est <- permTS(dep_values, nodep_values, alternative = "two.sided", 
                method = "exact.mc", control = permControl(p.conf.level = 0.95, nmc = 10000, seed = 123))
  
  # Store results
  results$p_value[i] <- est$p.value
  results$mean_diff[i] <- mean(dep_values) - mean(nodep_values)  # Median difference
  
  # Compute 95% confidence intervals (for whiskers)
  results$ci_lower_dep[i] <- quantile(dep_values, 0.025)
  results$ci_upper_dep[i] <- quantile(dep_values, 0.975)
  results$mean_dep[i] <- mean(dep_values)
  results$q1_dep[i] <- quantile(dep_values, 0.25)  # Q1 for boxplot
  results$q3_dep[i] <- quantile(dep_values, 0.75)  # Q3 for boxplot
  
  results$ci_lower_nodep[i] <- quantile(nodep_values, 0.025)
  results$ci_upper_nodep[i] <- quantile(nodep_values, 0.975)
  results$mean_nodep[i] <- mean(nodep_values)
  results$q1_nodep[i] <- quantile(nodep_values, 0.25)  # Q1 for boxplot
  results$q3_nodep[i] <- quantile(nodep_values, 0.75)  # Q3 for boxplot
}

# Apply multiple comparisons correction (Benjamini-Hochberg)
results$p_value_adj <- p.adjust(results$p_value, method = "BH")

# Format p-values for plotting
results$p_label <- ifelse(results$p_value_adj < 0.001, "<0.001", sprintf("p=%.3f", results$p_value_adj))

# Add median difference to p-value labels
results$p_label <- paste0(results$p_label, " (Δ-Median: ", round(results$mean_diff, 3), ")")

# Reshape data for boxplot visualization
plot_data <- data.frame(
  Predictor = rep(results$Predictor, each = 2),
  Mean = as.vector(t(cbind(results$mean_dep, results$mean_nodep))),
  Q1 = as.vector(t(cbind(results$q1_dep, results$q1_nodep))),
  Q3 = as.vector(t(cbind(results$q3_dep, results$q3_nodep))),
  Lower = as.vector(t(cbind(results$ci_lower_dep, results$ci_lower_nodep))),
  Upper = as.vector(t(cbind(results$ci_upper_dep, results$ci_upper_nodep))),
  Group = rep(c("Depressed", "Non-Depressed"), length(predictors))
)

library(dplyr)

# Create a new dataset with only the row with the highest Upper value for each Predictor
plot_data_single <- plot_data %>%
  group_by(Predictor) %>%
  filter(Upper == max(Upper)) %>%
  ungroup()

ggplot(plot_data, aes(x = Predictor, y = Mean, fill = Group)) +
  # Boxplot with proper dodge positioning
  geom_boxplot(aes(ymin = Lower, ymax = Upper, lower = Q1, upper = Q3, middle = Mean), 
               stat = "identity", position = position_dodge(width = 0.6), width = 0.5, alpha = 0.5) +
  
  geom_text(data = plot_data_single %>% filter(Predictor %in% c("CDQ009B", "CDQ009A", "CDQ009E")), 
            aes(
              x = Predictor, 
              y = Upper + 0.01,  # Offset above the upper value of the boxplot
              label = "*"),  # Label for significance
            size = 6, fontface = "bold", vjust = 0, family = "Times New Roman") +
  # Labels and theme adjustments
  labs(
    title = "Increased Predicted Probability of Coronary Artery Disease by Chest Pain Characterization\n Depressed vs Non-Depressed Cohorts\n
    *Sensitvity Analysis*",
    # subtitle = "Benjamini & Hochberg Corrected P-values: <0.001 ***",  # Subtitle with p-value
    x = "Predictor", 
    y = expression(Delta * "ALE / Increase in predicted probability of CAD (%)")
  ) +
  
  theme_minimal(base_family = "Times New Roman") +  # Change all text to Times New Roman
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, family = "Times New Roman"),  # Larger predictor labels
    axis.text.y = element_text(size = 12, family = "Times New Roman"),  
    axis.title.x = element_text(size = 12, family = "Times New Roman"),  
    axis.title.y = element_text(size = 12, family = "Times New Roman"),  
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold", family = "Times New Roman"),  # Center title
    plot.subtitle = element_text(size = 12, hjust = 0.5, family = "Times New Roman"),  # Center subtitle
    legend.text = element_text(size = 12, family = "Times New Roman"),
    legend.title = element_text(size = 12, family = "Times New Roman"),
    # Increase the left margin
    plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 30)  # Adjust the left margin to 30
  ) +
  
  scale_x_discrete(labels = new_labels) +  # Custom predictor labels
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))  # Convert y-axis to percent format


## Cohort specifc pain profile - within-group and across-group pairwise comparison

peripheral_dep <- rowMeans(cbind(diff_distributions_dep$CDQ009A, 
                                 diff_distributions_dep$CDQ009C,
                                 diff_distributions_dep$CDQ009E,
                                 diff_distributions_dep$CDQ009G,
                                 diff_distributions_dep$CDQ009H))

central_dep <- rowMeans(cbind(diff_distributions_dep$CDQ009D, diff_distributions_dep$CDQ009B,
                              diff_distributions_dep$CDQ009F
))


peripheral_nd<- rowMeans(cbind(diff_distributions_nodep$CDQ009A, diff_distributions_nodep$CDQ009C,
                               diff_distributions_nodep$CDQ009E,
                               diff_distributions_nodep$CDQ009G,
                               diff_distributions_nodep$CDQ009H))
central_nd<- rowMeans(cbind(diff_distributions_nodep$CDQ009D, diff_distributions_nodep$CDQ009B,
                            diff_distributions_nodep$CDQ009F)
)


diff_pct_dep <- quantile(peripheral_dep-central_dep, c(0.025,0.25,0.5,0.75, 0.975))
diff_pct_nd <- quantile(peripheral_nd-central_nd,c(0.025,0.25,0.5,0.75, 0.975))

diff_across_p <- quantile(peripheral_dep-peripheral_nd,c(0.025,0.25,0.5,0.75, 0.975))
diff_across_c <- quantile(central_dep-central_nd,c(0.025,0.25,0.5,0.75, 0.975))

comp <- data.frame()
comp <- rbind(diff_pct_dep, diff_pct_nd, diff_across_c, diff_across_p)
print(comp)
rownames(comp) <- c("Within-cohort: Depressed", "Within-cohort: Non-Depressed",
                    "Across-cohort: DCPP", "Across-cohort: NDCPP")


# Load necessary libraries
library(ggplot2)
library(reshape2)  # For reshaping data

# Convert matrix to long format for ggplot
comp_df <- as.data.frame(comp)
comp_df$Predictor <- rownames(comp_df)
comp_long <- melt(comp_df, id.vars = "Predictor")

# Boxplot with IQR, whiskers, and caps
ggplot(comp_long, aes(x = Predictor, y = value)) +
  geom_boxplot(fill = "lightblue", color = "black", width = 0.6, outlier.shape = NA) +  # Box with IQR, whiskers & caps
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +  # Dashed red line at y=0
  labs(title = "Within-cohort and Across-cohort differences in cohort-specifc pain profile ALE estimates\n *Sensitivity Analysis*", x = "Predictors", y = "Difference") +
  theme_minimal(base_family = "Times New Roman") +  # Set font to Times New Roman
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Center the title
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )
