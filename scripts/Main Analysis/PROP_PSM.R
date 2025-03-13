# Load required packages
library(survey)
library(tidyr)
library(ggplot2)

# List the eight binary outcome variables
vars <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D",
          "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")

### ----------------- 1. Compare Proportions by DEPR_NEW -----------------

# Calculate weighted proportions for the eight variables by DEPR_NEW
prop_depr <- svyby(formula = as.formula(paste("~", paste(vars, collapse = " + "))),
                   by = ~DEPR_NEW,
                   design = matched_design,
                   FUN = svymean,
                   keep.var = TRUE)
# --- Step 1: Build the summary data frame ---

# List the outcomes (base names without the level suffix)
outcomes <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", 
              "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")

# Initialize an empty list to store rows
plot_list <- list()

prop_depr_df <- as.data.frame(prop_depr)
# Loop over each outcome and each row of prop_depr_df
for(i in seq_along(outcomes)) {
  out <- outcomes[i]
  # For the "level 1" values, the mean column name is like "CDQ009A1"
  col_mean <- paste0(out, "1")
  # The corresponding SD (standard error) column is "se.CDQ009A1"
  col_sd <- paste0("se.", out, "1")
  
  # Loop over the two rows (one per DEPR_NEW group)
  for(j in 1:nrow(prop_depr_df)) {
    depr <- prop_depr_df$DEPR_NEW[j]
    mean_val <- prop_depr_df[j, col_mean]
    sd_val <- prop_depr_df[j, col_sd]
    # Calculate lower and upper bounds for the "box" (±1 SD) and "whiskers" (±2 SD)
    lower_box <- mean_val - sd_val
    upper_box <- mean_val + sd_val
    lower_whisker <- mean_val - 2 * sd_val
    upper_whisker <- mean_val + 2 * sd_val
    # Create a row with these values
    plot_list[[length(plot_list) + 1]] <- data.frame(
      Outcome = out,
      DEPR_NEW = as.factor(depr),
      Mean = mean_val,
      SD = sd_val,
      LowerBox = lower_box,
      UpperBox = upper_box,
      LowerWhisker = lower_whisker,
      UpperWhisker = upper_whisker
    )
  }
}

# Combine the list into one data frame
plot_df <- do.call(rbind, plot_list)

# --- Step 2: Plotting with ggplot2 ---
library(ggplot2)

# Update legend labels
plot_df$DEPR_NEW <- factor(plot_df$DEPR_NEW, levels = c(0, 1), labels = c("No", "Yes"))

# Use position_dodge() to place the DEPR_NEW groups side-by-side for each Outcome
p <- ggplot(plot_df, aes(x = Outcome, fill = DEPR_NEW)) +
  # Whiskers: thin error bars representing ±2 SD from the mean
  geom_errorbar(aes(ymin = LowerWhisker, ymax = UpperWhisker), 
                width = 0.2, 
                position = position_dodge(0.5)) +
  # "Box": thick errorbar representing ±1 SD from the mean 
  geom_errorbar(aes(ymin = LowerBox, ymax = UpperBox), 
                width = 0.5, 
                position = position_dodge(0.5), 
                size = 2) +
  # Mean as a point marker
  geom_point(aes(y = Mean), 
             position = position_dodge(0.5), 
             size = 4, shape = 21, color = "black") +  # Increased point size
  labs(x = "Outcome",
       y = "Proportion",
       fill = "Depression Status",
       title = "Chest Pain Location Survey-Weighted Proportions by Depression") +
  scale_x_discrete(labels = new_labels) +
  theme_minimal(base_family = "Times New Roman") +  # Set font to Times New Roman
  theme(
    text = element_text(family = "Times New Roman", size = 16),  # Set all text size
    axis.title = element_text(size = 18),  # Axis labels size
    axis.text = element_text(size = 16),  # Axis text size
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),  # Centered title
    legend.title = element_text(size = 18),  # Legend title size
    legend.text = element_text(size = 16),  # Legend text size
    legend.key.size = unit(2, "cm")  # Increase legend marker size
  ) 

print(p)


#### With chi square test

# --- Step 1: Build the summary data frame ---
# List the outcomes (base names without the level suffix)
outcomes <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", 
              "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")

# Initialize an empty list to store rows
plot_list <- list()

# Loop over each outcome and each row of prop_depr_df
for(i in seq_along(outcomes)) {
  out <- outcomes[i]
  # For the "level 1" values, the mean column name is like "CDQ009A1"
  col_mean <- paste0(out, "1")
  # The corresponding SD (standard error) column is "se.CDQ009A1"
  col_sd <- paste0("se.", out, "1")
  
  # Loop over the two rows (one per DEPR_NEW group)
  for(j in 1:nrow(prop_depr_df)) {
    depr <- prop_depr_df$DEPR_NEW[j]
    mean_val <- prop_depr_df[j, col_mean]
    sd_val <- prop_depr_df[j, col_sd]
    # Calculate lower and upper bounds for the "box" (±1 SD) and "whiskers" (±2 SD)
    lower_box <- mean_val - sd_val
    upper_box <- mean_val + sd_val
    lower_whisker <- mean_val - 2 * sd_val
    upper_whisker <- mean_val + 2 * sd_val
    # Create a row with these values
    plot_list[[length(plot_list) + 1]] <- data.frame(
      Outcome = out,
      DEPR_NEW = as.factor(depr),
      Mean = mean_val,
      SD = sd_val,
      LowerBox = lower_box,
      UpperBox = upper_box,
      LowerWhisker = lower_whisker,
      UpperWhisker = upper_whisker
    )
  }
}

# Combine the list into one data frame
plot_df <- do.call(rbind, plot_list)

# --- Step 2: Calculate chi-square test p-values for each outcome ---
# For each outcome, compare the binary variable (e.g. CDQ009A) across DEPR_NEW groups
# Note: Here we assume that in the survey design object, the variable is named as the outcome (e.g., "CDQ009A")
chi_pvals <- sapply(outcomes, function(var) {
  chi <- svychisq(as.formula(paste("~", var, "+ DEPR_NEW")), design = matched_design, statistic="adjWald")
  chi$p.value
})

# Build an annotation data frame with one row per outcome.
pvals_df <- data.frame(Outcome = outcomes, p_value = chi_pvals, stringsAsFactors = FALSE)

# For placing the p-value above the boxplots, find the maximum UpperWhisker for each outcome in plot_df.
annot_df <- aggregate(UpperWhisker ~ Outcome, data = plot_df, max)
# Merge with the p-values
annot_df <- merge(annot_df, pvals_df, by = "Outcome")
# Create a new column for the annotation y position (add a small offset)
annot_df$ypos <- annot_df$UpperWhisker + 0.02

# Format the p-values as a string (e.g., "p = 0.012")
annot_df$p_label <- paste0("p = ", formatC(annot_df$p_value, format = "f", digits = 3))

# --- Step 3: Plotting with ggplot2 ---
library(ggplot2)

plot_df$DEPR_NEW <- ifelse(plot_df$DEPR_NEW == 0, "Non-Depressed", "Depressed")

new_labels <- c(
  "CDQ009B" = "Right Chest Pain",
  "CDQ009D" = "Upper Sternal Pain",
  "CDQ009F" = "Left Chest Pain",
  "CDQ009H" = "Epigastric Pain",
  "CDQ009C" = "Neck Pain",
  "CDQ009E" = "Lower Sternal Pain",
  "CDQ009G" = "Left Arm Pain",
  "CDQ009A" = "Right Arm Pain"
)

library(ggplot2)

p <- ggplot(plot_df, aes(x = Outcome, fill = DEPR_NEW)) +
  # Whiskers: thin error bars representing ±2 SD from the mean
  geom_errorbar(aes(ymin = LowerWhisker, ymax = UpperWhisker), 
                width = 0.2, 
                position = position_dodge(0.5)) +
  # "Box": thick error bar representing ±1 SD from the mean 
  geom_errorbar(aes(ymin = LowerBox, ymax = UpperBox), 
                width = 0.5, 
                position = position_dodge(0.5), 
                size = 2) +
  # Mean as a point marker
  geom_point(aes(y = Mean), 
             position = position_dodge(0.5), 
             size = 3, shape = 21, color = "black") +
  # Add p-value text above each outcome
  geom_text(data = annot_df, 
            aes(x = Outcome, y = ypos, label = p_label, group = Outcome), 
            vjust = 0, size = 5, fontface = "italic", inherit.aes=FALSE) +
  # Update axis labels
  scale_x_discrete(labels = new_labels) +
  
  # Update fill scale and legend
  scale_fill_manual(values = c("Non-Depressed" = "lightblue", "Depressed" = "pink"), 
                    labels = c("Depressed", "Non-Depressed")) +
  
  # Update labels
  labs(x = "Pain Location",
       y = "Proportion",
       fill = "Group",  # Change legend title
       title = "Chest Pain Location Survey-Weighted Proportions \n and Adjusted-Wald Association by Depression") +
  
  theme_minimal(base_family = "Times New Roman") +  # Set global font to Times New Roman
  theme(
    text = element_text(family = "Times New Roman", size = 16),  # Increase overall text size
    axis.title.x = element_text(size = 18, face = "bold"),  # X-axis title size
    axis.title.y = element_text(size = 18, face = "bold"),  # Y-axis title size
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1),  # Rotate x-axis labels
    axis.text.y = element_text(size = 16),  # Y-axis labels size
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),  # Title size & centering
    legend.title = element_text(size = 18),  # Legend title size
    legend.text = element_text(size = 16)  # Legend text size
  )

print(p)





### ----------------- 2. Compare Proportions by DEPR_NEW, Stratified by CAD -----------------

# Calculate weighted proportions by both DEPR_NEW and CAD
prop_depr_cad <- svyby(formula = as.formula(paste("~", paste(vars, collapse = " + "))),
                       by = ~DEPR_NEW + CAD,
                       design = matched_design,
                       FUN = svymean,
                       keep.var = TRUE)

# Convert to data frame
prop_depr_cad_df <- as.data.frame(prop_depr_cad)

# Extract the estimates
prop_est_cad <- prop_depr_cad_df[, c("DEPR_NEW", "CAD", vars)]

# Extract standard errors
prop_se_cad <- prop_depr_cad_df[, grep("^se\\.", colnames(prop_depr_cad_df))]

# Clean the standard error column names
colnames(prop_se_cad) <- gsub("^se\\.", "", colnames(prop_se_cad))

# Convert to long format
prop_est_cad_long <- reshape(prop_est_cad, 
                             varying = vars, 
                             v.names = "prop", 
                             times = vars, 
                             timevar = "variable", 
                             direction = "long")

prop_se_cad_long <- reshape(prop_se_cad, 
                            varying = colnames(prop_se_cad), 
                            v.names = "se", 
                            times = colnames(prop_se_cad), 
                            timevar = "variable", 
                            direction = "long")

# Merge the estimates and standard errors
prop_final_cad <- merge(prop_est_cad_long, prop_se_cad_long, by = c("DEPR_NEW", "CAD", "variable"))

# Plot: Bar plot with error bars stratified by CAD, faceted by DEPR_NEW
p2 <- ggplot(prop_final_cad, aes(x = variable, y = prop, fill = factor(CAD))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = prop - se, ymax = prop + se),
                width = 0.2,
                position = position_dodge(width = 0.9)) +
  facet_wrap(~ DEPR_NEW) +
  labs(x = "Variable",
       y = "Weighted Proportion",
       fill = "CAD") +
  theme_minimal() +
  ggtitle("Weighted Proportions Stratified by CAD within DEPR_NEW Groups")

print(p2)




##### Number of chest pain location by depression

library(survey)
library(ggplot2)
library(ggpubr)
library(extrafont)

# Convert binary factor variables to numeric and sum them
matched_design <- update(matched_design, 
                         total_pain_regions = rowSums(sapply(matched_design$variables[, c("CDQ009A", "CDQ009B", "CDQ009C", 
                                                                                          "CDQ009D", "CDQ009E", "CDQ009F", 
                                                                                          "CDQ009G", "CDQ009H")], 
                                                             function(x) as.numeric(as.character(x))), na.rm = TRUE))
update_design <- subset(matched_design, matched_design$variables$total_pain_regions > 0)

# Compute survey-weighted mean and 95% confidence intervals
svy_mean <- svyby(~total_pain_regions, ~DEPR_NEW, update_design, svymean, vartype = c("ci"))
colnames(svy_mean) <- c("DEPR_NEW", "mean", "ci_l", "ci_u")  # Rename columns for clarity

# Perform a weighted t-test
svy_ttest <- svyttest(total_pain_regions ~ DEPR_NEW, update_design)

# Extract formatted p-value
p_value <- formatC(svy_ttest$p.value, digits = 3)

update_df <- update_design$variables

# Create boxplot with weighted means and confidence intervals
ggplot(update_df, aes(x = as.factor(DEPR_NEW), y = total_pain_regions)) +
  geom_point(data = svy_mean, aes(x = as.factor(DEPR_NEW), y = mean), 
             color = "red", size = 4) +  # Weighted means
  geom_errorbar(data = svy_mean, aes(x = as.factor(DEPR_NEW), ymin = ci_l, ymax = ci_u), 
                width = 0.2, color = "darkblue", size = 1.2, inherit.aes = FALSE) +  # Thicker, dark blue error bars
  labs(title = "Survey-weighted Total Chest Pain Regions by Depression Status",
       x = "Depression Status", 
       y = "Weighted Mean Total Pain Regions") +
  annotate("text", x = 1.5, y = 1.825, family = "Times New Roman",
           label = paste("p =", p_value), size = 8, color = "black") +
  theme_minimal() +
  scale_x_discrete(labels = c("No Depression", "Depression")) +
  theme(
    text = element_text(family = "Times New Roman", size = 14, hjust=0.5),  # Set font to Times New Roman and increase text size
    plot.title = element_text(size = 18, face = "bold", hjust=0.5),  # Increase title size and make it bold
    axis.title.x = element_text(size = 16),  # Increase x-axis label size
    axis.title.y = element_text(size = 16),  # Increase y-axis label size
    axis.text.x = element_text(size = 14),  # Increase x-axis tick label size
    axis.text.y = element_text(size = 14)   # Increase y-axis tick label size
  )



#### Survey weighted correlation

library(survey)
library(ggplot2)
library(reshape2)

# Define pain region variables (adjust if needed)
pain_vars <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", 
               "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")

# Convert CDQ009A-CDQ009H to numeric (0/1)
pain_vars <- paste0("CDQ009", LETTERS[1:8])  # Adjust if needed


new_design <- matched_design

# Convert variables in the survey design object
new_design$variables[pain_vars] <- lapply(
  new_design$variables[pain_vars], 
  function(x) as.numeric(as.character(x)))  # Handles factors safely

# Subset designs by depression status
depr0_design <- subset(new_design, DEPR_NEW == 0)
depr1_design <- subset(new_design, DEPR_NEW == 1)

# Function to compute correlation matrix for a survey design
compute_svy_cor <- function(design, variables) {
  # Compute covariance matrix
  cov_mat <- svyvar(as.formula(paste("~", paste(variables, collapse = "+"))), 
                    design = design, na.rm = TRUE)
  
  # Extract covariance matrix
  cov_mat <- as.matrix(cov_mat)
  
  # Compute correlation matrix
  cor_mat <- cov2cor(cov_mat)
  
  return(cor_mat)
}

# Compute correlation matrices for both groups
cor_depr0 <- compute_svy_cor(depr0_design, pain_vars)
cor_depr1 <- compute_svy_cor(depr1_design, pain_vars)

# Melt correlation matrices into data frames
cor_df0 <- melt(cor_depr0)
cor_df0$group <- "No Depression"
cor_df1 <- melt(cor_depr1)
cor_df1$group <- "Depression"

# Combine data
combined_cor <- rbind(cor_df0, cor_df1)
colnames(combined_cor) <- c("Var1", "Var2", "Phi", "Group")

ggplot(combined_cor, aes(x = Var1, y = Var2, fill = Phi)) +
  geom_tile(color = "white") +
  facet_wrap(~Group) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", 
    midpoint = 0, limits = c(-1, 1),
    name = "Phi Coefficient"
  ) +
  labs(
    title = "Survey-Weighted Phi Correlation Coefficients Between Pain Regions",
    subtitle = "Stratified by Depression Status",
    x = "", 
    y = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(family = "Times New Roman", size = 16),
    plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 16, hjust=0.5)
  ) +
  scale_x_discrete(labels=new_labels) +
  scale_y_discrete(labels=new_labels)




### Correlation test

library(survey)
library(ggplot2)
library(dplyr)
library(tidyr)

# Convert binary variables to 0/1 (if needed)
pain_vars <- paste0("CDQ009", LETTERS[1:8])  # CDQ009A to CDQ009H
new_design$variables[pain_vars] <- lapply(
  new_design$variables[pain_vars], 
  function(x) as.numeric(as.character(x))
)

# Create replicate design
rep_design <- as.svrepdesign(new_design, type = "bootstrap", replicates = 1000)

# Function to compute weighted Phi for a single replicate
compute_phi_replicate <- function(w, data, var1, var2) {
  data$weights <- w
  mean1 <- weighted.mean(data[[var1]], data$weights, na.rm = TRUE)
  mean2 <- weighted.mean(data[[var2]], data$weights, na.rm = TRUE)
  cov <- sum(data$weights * (data[[var1]] - mean1) * (data[[var2]] - mean2)) / sum(data$weights)
  var1 <- sum(data$weights * (data[[var1]] - mean1)^2) / sum(data$weights)
  var2 <- sum(data$weights * (data[[var2]] - mean2)^2) / sum(data$weights)
  phi <- cov / sqrt(var1 * var2)
  return(phi)
}

# Generate all unique pairs of pain regions
pain_vars <- paste0("CDQ009", LETTERS[1:8])
pairs <- combn(pain_vars, 2, simplify = FALSE)

# Initialize empty list to store results
results <- list()

# Loop through all pairs
for (pair in pairs) {
  var1 <- pair[1]
  var2 <- pair[2]
  
  # Compute Phi coefficients for both groups
  phi0 <- withReplicates(depr0_rep, compute_phi_replicate, var1 = var1, var2 = var2)
  phi1 <- withReplicates(depr1_rep, compute_phi_replicate, var1 = var1, var2 = var2)
  
  # Extract estimates and SEs
  phi0_est <- phi0[1]
  phi1_est <- phi1[1]
  se0 <- sqrt(attr(phi0, "var"))
  se1 <- sqrt(attr(phi1, "var"))
  
  # Compute difference statistics
  diff_phi <- phi1_est - phi0_est
  se_diff <- sqrt(se0^2 + se1^2)
  z_score <- diff_phi / se_diff
  p_value <- 2 * pnorm(-abs(z_score))
  
  # Store results
  results[[paste(var1, var2, sep = "-")]] <- data.frame(
    Var1 = var1,
    Var2 = var2,
    Phi_NoDepression = phi0_est,
    Phi_Depression = phi1_est,
    Difference = diff_phi,
    SE = se_diff,
    Z = z_score,
    p_value = p_value
  )
}

# Combine all results into a dataframe
results_df <- bind_rows(results)

# Bonferroni correction
bf_sig <- 0.05/28


# Create symmetric matrix for plotting
plot_df <- results_df %>%
  dplyr::select(Var1, Var2, Difference, p_value) %>%
  dplyr::mutate(
    p_text = case_when(
      p_value < bf_sig ~ "***",
      TRUE ~ format(round(p_value, 2), nsmall = 2)
    )
  )

# Create heatmap with p-value annotations
ggplot(plot_df, aes(x = Var1, y = Var2, fill = Difference)) +
  geom_tile(color = "white") +
  geom_text(aes(label = p_text), color = "black", size = 4) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0, limits = c(-1, 1),
    name = expression(Delta*"Phi")
  ) +
  labs(
    title = "Difference in Phi Correlation Coefficients Between Depression Cohorts",
    subtitle = "Bonferroni Corrected Significance Level: ~0.0018",
    x = "",
    y = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(family = "Times New Roman", size = 18),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5)
  ) +
  coord_fixed() + # Keep aspect ratio square
  scale_x_discrete(labels=new_labels) +
  scale_y_discrete(labels=new_labels)