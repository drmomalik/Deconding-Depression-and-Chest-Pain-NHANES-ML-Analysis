## How to use survery package to do analysis with weights from NHANES
library(survey)
library(sf)
library(ggplot2)
library(dplyr)
library(ggpattern)

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
model <- svyglm(CDQ009A ~ DEPR_LVL, 
                design = nhanes_design, 
                family = quasibinomial())  # Use quasibinomial for survey-adjusted logit
summary(model)


## Log regression for all levels of depression and body regions 

# Template data
body_parts <- c("right_arm", "right_chest", "neck", "upper_sternum", 
                "lower_sternum", "left_chest", "left_arm", "epigastric")
variables <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", 
               "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")


# Placeholder for results
hbdata_list <- list()

# Iterate over depression levels (1 to 5)
for (depr_lvl in 1:5) {
  # Create a new dataframe for this depression level
  hbdata_new <- hbdata
  
  # Iterate over body region variables
  for (i in seq_along(variables)) {
    # Fit logistic regression for the current body region
    formula <- as.formula(paste0(variables[i], " ~ -1 + DEPR_LVL"))
    model <- svyglm(formula, design = nhanes_design, family = quasibinomial())
    
    # Extract coefficients and p-values
    coef_vals <- coef(model)
    p_vals <- summary(model)$coefficients[, "Pr(>|t|)"]
    
    # Calculate odds ratio for the current depression level
    odds_ratio <- exp(coef_vals[paste0("DEPR_LVL", depr_lvl)])
    
    # Determine significance for the current coefficient
    significant <- ifelse(p_vals[paste0("DEPR_LVL", depr_lvl)] < 0.05, 1, 0)
    
    # Update the hbdata template
    hbdata_new$value[hbdata_new$body_part == body_parts[i]] <- odds_ratio
    hbdata_new$significant[hbdata_new$body_part == body_parts[i]] <- significant
  }
  
  # Save the updated dataframe
  hbdata_list[[depr_lvl]] <- hbdata_new
}

# Combine all dataframes for plotting
human_body_data <- bind_rows(
  hbdata_list %>%
    purrr::imap(~ .x %>% mutate(DEPR_LVL = .y))
)

# Plot faceted heatmaps with significance-based alpha and bold outlines
ggplot(human_body_data) +
  geom_sf(aes(fill = value, alpha = significant), color = "black", linewidth = ifelse(human_body_data$significant == 1, 0.8, 0.2)) +
  scale_fill_gradient(low = "blue", high = "red") +
  scale_alpha_continuous(range = c(0.3, 1)) +  # Fade non-significant regions
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  ) +
  labs(fill = "Proportion", alpha = "Significance") +
  facet_wrap(~DEPR_LVL)


## Alt: Plot using DEPR_LVL1 as the reference/baseline 

hbdata <- human_body_data

# Template data
body_parts <- c("right_arm", "right_chest", "neck", "upper_sternum", 
                "lower_sternum", "left_chest", "left_arm", "epigastric")
variables <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", 
               "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")


# Placeholder for results
hbdata_list <- list()

# Iterate over depression levels (2 to 5)
for (depr_lvl in 2:5) {
  # Create a new dataframe for this depression level
  hbdata_new <- hbdata
  
  # Iterate over body region variables
  for (i in seq_along(variables)) {
    # Fit logistic regression for the current body region
    formula <- as.formula(paste0(variables[i], " ~ DEPR_LVL"))
    model <- svyglm(formula, design = nhanes_design, family = quasibinomial())
    
    # Extract coefficients and p-values
    coef_vals <- coef(model)
    p_vals <- summary(model)$coefficients[, "Pr(>|t|)"]
    
    # Calculate odds ratio for the current depression level
    odds_ratio <- exp(coef_vals[paste0("DEPR_LVL", depr_lvl)])
    
    # Determine significance for the current coefficient
    significant <- ifelse(p_vals[paste0("DEPR_LVL", depr_lvl)] < 0.05, 1, 0)
    
    # Update the hbdata template
    hbdata_new$value[hbdata_new$body_part == body_parts[i]] <- odds_ratio
    hbdata_new$significant[hbdata_new$body_part == body_parts[i]] <- significant
  }
  
  # Save the updated dataframe
  hbdata_list[[depr_lvl]] <- hbdata_new
}

# Assign dataframes to individual variables for each depression level
for (depr_lvl in 1:5) {
  assign(paste0("hbdata_", depr_lvl), hbdata_list[[depr_lvl]])
}

# Add DEPR_LVL column to each dataframe and combine using rbind
hbdata_2 <- hbdata_2 %>% mutate(DEPR_LVL = "2")
hbdata_3 <- hbdata_3 %>% mutate(DEPR_LVL = "3")
hbdata_4 <- hbdata_4 %>% mutate(DEPR_LVL = "4")
hbdata_5 <- hbdata_5 %>% mutate(DEPR_LVL = "5")

# Combine all dataframes
human_body_data <- rbind(hbdata_2, hbdata_3, hbdata_4, hbdata_5)

# Add custom attributes for plotting
human_body_data <- human_body_data %>%
  mutate(
    pattern = ifelse(significant == 1, "stripe", "none"),  # Pattern for significant
    alpha = ifelse(significant == 1, 1, 0.2)                  # Transparency for non-significant
  )

# Plot faceted heatmaps
ggplot(human_body_data) +
  geom_sf(aes(fill = value), color = "black", linewidth = 1) + 
  geom_sf_pattern(aes(fill = value, pattern = pattern, alpha = alpha), color = "black",
                 pattern_density = 0.2) +
  scale_fill_gradient(low = "blue", high = "red") +
  scale_alpha_continuous(range = c(0.2, 1), guide = "none") +  # Fade non-significant regions
  scale_pattern_manual(values = c("none", "stripe")) +     # Specify patterns
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  ) +
  labs(fill = "Odds-Ratio", pattern = "Significance") +
  facet_wrap(~DEPR_LVL)




### Alt: Depr as binary (=>15) and body region


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

# subset design object
subset_data <- subset(
  nhanes_design,  # Access the data (variables) in the svydesign object
  !apply(
    nhanes_design$variables[, c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")],
    1, 
    function(x) all(x == 0 | is.na(x))  # Keep rows where not all values of CDQ009A-CDQ009H are 0 or NA
  ))

hbdata <- human_body_data

# Template data
body_parts <- c("right_arm", "right_chest", "neck", "upper_sternum", 
                "lower_sternum", "left_chest", "left_arm", "epigastric")
variables <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", 
               "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")


# Placeholder for results
hbdata_new <- hbdata

# Iterate over body region variables
  for (i in seq_along(variables)) {
    # Fit logistic regression for the current body region
    formula <- as.formula(paste0(variables[i], " ~ -1+DEPR_BIN"))
    model <- svyglm(formula, design = subset_data, family = quasibinomial())
    
    # Extract coefficients and p-values
    coef_vals <- coef(model)
    p_vals <- summary(model)$coefficients[, "Pr(>|t|)"]
    
    # Calculate odds ratio for the current depression level
    odds_ratio <- exp(coef_vals[paste0("DEPR_BIN", depr_lvl)])
    
    # Determine significance for the current coefficient
    significant <- ifelse(p_vals[paste0("DEPR_BIN", depr_lvl)] < 0.05, 1, 0)
    
    # Update the hbdata template
    hbdata_new$value[hbdata_new$body_part == body_parts[i]] <- odds_ratio
    hbdata_new$significant[hbdata_new$body_part == body_parts[i]] <- significant
  }


# Combine all dataframes
hbd <- hbdata_new

# Add custom attributes for plotting
hbd <- hbd %>%
  mutate(
    pattern = ifelse(significant == 1, "stripe", "none"),  # Pattern for significant
    alpha = ifelse(significant == 1, 1, 0.2)                  # Transparency for non-significant
  )

# Plot faceted heatmaps
ggplot(hbd) +
  geom_sf(aes(fill = value), color = "black", linewidth = 1) + 
  geom_sf_pattern(aes(fill = value, pattern = pattern, alpha = alpha), color = "black",
                  pattern_density = 0.2) +
  scale_fill_gradient(low = "blue", high = "red") +
  scale_alpha_continuous(range = c(0.2, 1), guide = "none") +  # Fade non-significant regions
  scale_pattern_manual(values = c("none", "stripe")) +     # Specify patterns
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  ) +
  labs(fill = "Odds-Ratio", pattern = "Significance")