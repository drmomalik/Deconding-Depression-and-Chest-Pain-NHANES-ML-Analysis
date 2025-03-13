library(tidyverse)
library(nhanesA)
library(dplyr)


## How many observations of chest pain by location
loc_total <- numeric(8)
for (i in 18:25) {
  loc_total[i-17] <- length(df_lim$SEQN) - sum(is.na(df_lim[i]))
}
print(loc_total)


## How many respondents had chest pain in multiple locations 
df_lim$tot_loc_cp <- rowSums(df_lim[, c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", 
                                        "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")])
counts <- table(df_lim$tot_loc_cp)
print(counts)

colnames(allyr_df_coded)[8:15] <- c("Right Arm", "Right Chest",
                          "Neck", "Upper Sternum",
                          "Lower Sternum", "Left Chest",
                          "Left arm", "Epigastric")

## Sample histograms of depression level and chest pain location
library(ggplot2)

# Reshape the data to long format for CDQ009A through CDQ009H
allyr_long <- allyr_df_coded %>%
  pivot_longer(cols = 8:15, names_to = "Outcome", values_to = "BinaryOutcome")

# Summarize data to get counts and proportions for each level of DEPR_LVL and each Outcome
allyr_summary <- allyr_long %>%
  group_by(DEPR_LVL, Outcome, BinaryOutcome) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(DEPR_LVL, Outcome) %>%
  mutate(total_count = sum(count),
         proportion = count / total_count * 100) %>%
  ungroup()

# Plot the histograms with faceting by Outcome
ggplot(allyr_summary, aes(x = DEPR_LVL, y = count, fill = factor(BinaryOutcome))) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0("Total: ", total_count), 
                y = total_count + 1), 
            vjust = -0.5, size = 3, color = "black", 
            data = distinct(allyr_summary, DEPR_LVL, Outcome, .keep_all = TRUE)) +
  geom_text(aes(label = paste0(round(proportion, 1), "%")), 
            position = position_stack(vjust = 0.5), size = 3, color = "white") +
  labs(x = "DEPR_LVL", y = "Count", fill = "Binary Outcome") +
  facet_wrap(~ Outcome) +
  theme_minimal()


## To plot only proportions


# Reshape the data to long format for CDQ009A through CDQ009H
allyr_long <- allyr_df_coded %>%
  pivot_longer(cols = 8:15, names_to = "Outcome", values_to = "BinaryOutcome")

# Summarize data to get counts and proportions for each level of DEPR_LVL and each Outcome (only for 1's)
allyr_summary <- allyr_long %>%
  group_by(DEPR_LVL, Outcome, BinaryOutcome) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(DEPR_LVL, Outcome) %>%
  mutate(total_count = sum(count),
         proportion_ones = ifelse(BinaryOutcome == 1, count / total_count * 100, 0)) %>%
  ungroup()

# Calculate min and max proportions for each Outcome
proportion_range <- allyr_summary %>%
  group_by(Outcome) %>%
  summarise(min_prop = min(proportion_ones[proportion_ones > 0], na.rm = TRUE),
            max_prop = max(proportion_ones, na.rm = TRUE), .groups = "drop")

# Merge proportion_range with the main summary data
allyr_summary <- allyr_summary %>%
  left_join(proportion_range, by = "Outcome") %>%
  mutate(
    normalized_fill = (proportion_ones - min_prop) / (max_prop - min_prop) # Normalize proportions
  )

# Plot with individual scales per facet
ggplot(allyr_summary %>% filter(BinaryOutcome == 1), 
       aes(x = DEPR_LVL, y = proportion_ones, fill = normalized_fill)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(round(proportion_ones, 1), "%"), 
                y = proportion_ones / 2), # Adjust label position to be centered in the bar
            size = 3, color = "white") +
  labs(x = "DEPR_LVL", y = "Proportion of 1's (%)") +
  facet_wrap(~ Outcome, scales = "free_y") + 
  scale_y_continuous(labels = scales::percent_format(scale = 1)) + # Format as percentage
  scale_fill_gradient(low = "blue", high = "red", name = "Proportion (%)") +
  theme_minimal() +
  theme(strip.text = element_text(size = 10), # Adjust facet label size
        panel.spacing = unit(1, "limes"), # Increase space between panels
        axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels if needed


## Unadjusted simple log models for PHQ-9 score and pain location
df <- allyr_df %>%
  mutate(across(c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", 
                  "CDQ009F", "CDQ009G", "CDQ009H", "CDQ008", "CDQ010"), 
                ~ ifelse(is.na(.), 0, 1)))
# Create a vector of outcome variable names
outcomes <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", 
              "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")

# Initialize an empty dataframe to store results
odds_ratios_df <- data.frame()

# Loop through each outcome variable, fit the model, and calculate odds ratios
for (outcome in outcomes) {
  # Fit the model
  model <- glm(as.formula(paste(outcome, "~ -1 + DEPR_LVL")), family = binomial, data = df)
  
  # Extract coefficients and calculate odds ratios
  coef_data <- data.frame(
    DEPR_LVL = names(coef(model)),          # Levels of DEPR_LVL
    Odds_Ratio = exp(coef(model)),          # Convert log-odds to odds ratios
    Outcome = outcome                       # Store the outcome name
  )
  
  # Combine with the main dataframe
  odds_ratios_df <- rbind(odds_ratios_df, coef_data)
}

# Clean DEPR_LVL to remove the prefix (if it exists in your data)
odds_ratios_df$DEPR_LVL <- gsub("DEPR_LVL", "", odds_ratios_df$DEPR_LVL)

# Define custom facet titles
facet_titles <- c(
  CDQ009A = "Right Arm",
  CDQ009B = "Right Chest",
  CDQ009C = "Neck",
  CDQ009D = "Upper Sternum",
  CDQ009E = "Lower Sternum",
  CDQ009F = "Left Chest",
  CDQ009G = "Left Arm",
  CDQ009H = "Epigastric"
)

# Plot the odds ratios with custom facet titles
ggplot(odds_ratios_df, aes(x = DEPR_LVL, y = Odds_Ratio)) +
  geom_point(size = 3, color = "blue") +
  geom_lime(group = 1, color = "blue") +
  scale_y_log10() + # Log scale for odds ratios
  labs(
    x = "DEPR_LVL",
    y = "Odds Ratio",
    title = "Odds Ratios by Depression Level"
  ) +
  facet_wrap(~ Outcome, scales = "free_y", labeller = as_labeller(facet_titles)) + # Apply custom titles
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),  # Adjust facet label size
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    panel.spacing = unit(1, "limes")  # Increase space between panels
  )

## Examine trends in location of pain and depression lvl across years

# Load the required packages
library(dplyr)

# Specify the columns of interest
columns_of_interest <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", 
                         "CDQ009F", "CDQ009G", "CDQ009H", "DEPR_LVL")

# Create a data frame to store the chi-square test results
chi_results_df <- data.frame(Variable = character(),
                             Chi_Square_Stat = numeric(),
                             p_value = numeric(),
                             stringsAsFactors = FALSE)

# Loop through each variable and perform a chi-square test
for (var in columns_of_interest) {
  
  # Create a contingency table for the variable against SDDSRVYR
  contingency_table <- table(df_lim[[var]], df_lim$SDDSRVYR)
  
  # Perform the chi-square test
  chi_test <- chisq.test(contingency_table)
  
  # Store the results
  chi_results_df <- chi_results_df %>%
    add_row(Variable = var,
            Chi_Square_Stat = chi_test$statistic,
            p_value = chi_test$p.value)
}

# View the chi-square test results
print(chi_results_df)

## Plot proportions of variable of interest over cycles with chi-sqaure

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Specify the columns of interest
columns_of_interest <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", 
                         "CDQ009F", "CDQ009G", "CDQ009H", "DEPR_LVL")

# Create a vector of descriptive titles
variable_titles <- c(
  "CDQ009A" = "Right Arm Pain",
  "CDQ009B" = "Right Chest Pain",
  "CDQ009C" = "Neck Pain",
  "CDQ009D" = "Upper Sternum Pain",
  "CDQ009E" = "Lower Sternum Pain",
  "CDQ009F" = "Left Chest Pain",
  "CDQ009G" = "Left Arm Pain",
  "CDQ009H" = "Epigastric Pain",
  "DEPR_LVL" = "Depression Level Pain"
)

# Function to plot proportions for each variable
plot_proportions <- function(variable, df) {
  
  # Calculate proportions for each level of the variable by SDDSRVYR
  df_long <- df %>%
    group_by(SDDSRVYR, !!sym(variable)) %>%
    tally() %>%
    group_by(SDDSRVYR) %>%
    mutate(proportion = n / sum(n)) %>%
    ungroup()
  
  # Create contingency table
  contingency_table <- table(df[[variable]], df$SDDSRVYR)
  
  # Perform chi-square test
  chi_test <- chisq.test(contingency_table)
  
  # If expected frequencies are too small, use Fisher's Exact Test instead
  if (any(chi_test$expected < 5)) {
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    p_value <- fisher_test$p.value
    chi_stat <- NULL  # Fisher's test does not return a chi-square statistic
  } else {
    p_value <- chi_test$p.value
    chi_stat <- chi_test$statistic  # Chi-square statistic for chi-square test
  }
  
  # Format p-value
  p_value_formatted <- ifelse(p_value < 0.01, "<0.01", round(p_value, 3))
  
  # Add chi-square/Fisher's results as text annotation
  if (!is.null(chi_stat)) {
    chi_text <- paste0("Chi^2: ", round(chi_stat, 2), "\nP: ", p_value_formatted)
  } else {
    chi_text <- paste0("P: ", p_value_formatted)  # Only show p-value for Fisher's test
  }
  
  # Get the descriptive title for the variable
  plot_title <- variable_titles[variable]
  
  # Create the line plot
  ggplot(df_long, aes(x = as.factor(SDDSRVYR), y = proportion, color = as.factor(!!sym(variable)))) +
    geom_line(aes(group = !!sym(variable)), size = 1) +
    geom_point() +
    labs(title = paste("Proportions of", plot_title, "by Cycle"),
         x = "Cycle", y = "Proportion", color = "Levels of Variable") +
    annotate("text", x = 1, y = max(df_long$proportion, na.rm = TRUE), 
             label = chi_text, hjust = 0, size = 4) +
    theme_minimal()
}

# Loop through each variable and create the plots
for (var in columns_of_interest) {
  print(plot_proportions(var, df_lim))
}

# Plot proportion of chest pain vs cycle
var <- "CDQ001"
variable_titles <- c("CDQ001" = "Chest Pain")
print(plot_proportions(var, df))





