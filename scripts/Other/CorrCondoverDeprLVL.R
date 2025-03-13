library(survey)
library(corrplot)
library(ggplot2)
library(dplyr)
library(jtools)


# List of binary variables
binary_vars <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")

# Define the survey design
sub_design <- svydesign(id = ~SDMVPSU, weights = ~MEC15YR,
                        strata = ~SDMVSTRA, data = dfc, nest = TRUE)

# Subset the survey design based on conditions for DEPR_LVL
sub_design <- subset(sub_design, !apply(
  sub_design$variables[, binary_vars],
  1, 
  function(x) all(x == 0 | is.na(x))  # Keep rows where not all values of CDQ009A-CDQ009H are 0 or NA
))

# Initialize a list to store the results
correlation_results <- list()

# Loop through each DEPR_LVL (1 to 5)
for (depr_lvl in 1:5) {
  
  # Subset the design based on DEPR_LVL
  sub_design_lvl <- subset(sub_design, DEPR_LVL == depr_lvl)
  
  # Calculate the correlation matrix for the current depression level
  corr_result <- svycor(
    ~CDQ009A + CDQ009B + CDQ009C + CDQ009D + CDQ009E + CDQ009F + CDQ009G + CDQ009H,
    design = sub_design_lvl,
    na.rm = TRUE,
    sig.stats = TRUE,
    bootn = 500
  )
  
  # Extract the correlation matrix and p-values
  corr_matrix <- corr_result$cor
  p_matrix <- corr_result$p.values
  
  # Check if the correlation matrix is 8x8
  if (nrow(corr_matrix) == 8 && ncol(corr_matrix) == 8) {
    # Loop through each pair of body regions and store the correlation values and p-values
    for (i in 1:8) {
      for (j in setdiff(1:8, i)) {  # Only consider the upper triangle, avoiding duplicate pairs
        correlation_results[[length(correlation_results) + 1]] <- data.frame(
          DEPR_LVL = depr_lvl,
          Body_Region_1 = colnames(corr_matrix)[i],
          Body_Region_2 = colnames(corr_matrix)[j],
          Correlation = corr_matrix[i, j],
          p_value = p_matrix[i, j]
        )
      }
    }
  } else {
    warning(paste("Correlation matrix for DEPR_LVL", depr_lvl, "is not 8x8."))
  }
}

# Convert correlation results to a data frame for ggplot
plot_data <- do.call(rbind, correlation_results) %>%
  mutate(
    body_region = factor(Body_Region_1, levels = names(body_region_labels), labels = body_region_labels),
    other_region = factor(Body_Region_2, levels = names(body_region_labels), labels = body_region_labels),
    significant = p_value < 0.05  # Flag significant correlations
  )

# Plot the results
ggplot(plot_data, aes(x = DEPR_LVL, y = Correlation, color = other_region, group = other_region)) +
  geom_point(aes(shape = significant), size = 3) +  # Points with shape for significance
  scale_shape_manual(values = c(16, 8)) +  # Circle for non-significant, star for significant
  geom_smooth(
    aes(group = other_region), 
    method = "loess", se = TRUE, size = 1
  ) +  # Smoothed lines
  scale_color_manual(values = RColorBrewer::brewer.pal(8, "Set2")) +
  facet_wrap(~body_region, nrow = 2, ncol = 4, scales = "free_y") +
  labs(
    title = "Correlation between Chest Pain Regions Across Depression Levels",
    x = "Depression Level",
    y = "Non-parametric Bootstrap\nPhi coefficient",
    color = "Chest Pain Regions",
    shape = "Significant"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  )
