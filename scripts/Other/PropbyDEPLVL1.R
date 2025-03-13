library(survey)
library(dplyr)
library(purrr)
library(ggplot2)
library(sf)
library(progress)
library(mice)
library(future)
library(future.apply)

# Set up a parallel backend (e.g., 4 cores)
plan(multisession, workers = 4)  

# Subset variables to improve computational efficiency
responses <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")
predictors <- c("DEPR_TOT", "DEPR_BIN", "DEPR_LVL", "AGE_BIN", "RIAGENDR", "RIDRETH1", "SMQ020",
                "SMQ040", "HIQ011", "BPQ020", "BMXBMI", "BMI_LVL", "BPQ080",
                "ALQ130", "MEDDEP", "DMDBORNT", "PAQMV", "CADTOT", "DIDTOT",
                "DUQTOT", "INC_BIN", "INC3", "DMDEDUC2", "CDQ008", "CDQ010", "TOT_REG")
wd_subset <- wd[, c(responses, predictors,"SDDSRVYR", "SDMVPSU", "SDMVSTRA", "MEC15YR", "SEQN")]

# Initialize the predictor matrix for MICE
predictorMatrix <- mice::make.predictorMatrix(wd_subset)
predictorMatrix[, !colnames(predictorMatrix) %in% predictors] <- 0
predictorMatrix[!rownames(predictorMatrix) %in% predictors, ] <- 0

# Run MICE
m <- 5
maxit <- 5
progress_bar <- progress_bar$new(total = m*maxit)
imputations <- mice(wd_subset, 
                    seed = 123, 
                    maxit = maxit, 
                    m = m, 
                    predictorMatrix = predictorMatrix, 
                    printFlag = TRUE,
                    parallel = "future")
# Reset to sequential processing when done
plan(sequential)


## Run loop to plot proportions in heat map

# Template data
body_parts <- c("right_arm", "right_chest", "neck", "upper_sternum", 
                "lower_sternum", "left_chest", "left_arm", "epigastric")
variables <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", 
               "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")


# Assign template shapefile 
hbdata <- hbd

# Placeholder for results
hbdata_list <- list()

# Iterate over depression levels (1 to 5)
for (depr_lvl in 1:5) {
  # Create a new dataframe for this depression level
  hbdata_new <- hbdata
  
  # Iterate over body region variables (CDQ009A to CDQ009H)
  for (i in seq_along(variables)) {
    # Placeholder for proportions and variances across imputations
    prop_values <- numeric()
    prop_variances <- numeric()
    
    # Process each imputed dataset
    dat_mice <- complete(imputations, "all")
    for (dat in dat_mice) {
      # Define the survey design for this dataset
      nhanes_design <- svydesign(
        id = ~SDMVPSU,
        weights = ~MEC15YR,
        strata = ~SDMVSTRA,
        data = dat,
        nest = TRUE
      )
      #Subset for those who responded to questionnaire of CP
      svy_subset1 <- subset(nhanes_design, !apply(
        nhanes_design$variables[, responses], 1,
        function(x) all(x == 0 | is.na(x))
      ))
      
      # Subset the survey data to the current depression level
      svy_subset <- subset(svy_subset1, DEPR_LVL == depr_lvl)
      
      # Calculate survey-weighted proportion for the current body region variable (e.g., CDQ009A, CDQ009B)
      prop <- svymean(as.formula(paste0("~", variables[i])), svy_subset, na.rm = TRUE)
      
      # Extract the proportion and standard error for '1' (second row)
      prop_values <- c(prop_values, as.numeric(prop[2]))  # prop[2] gives the proportion for '1'
      prop_variances <- c(prop_variances, attr(prop, "var")[2])  # Standard error for '1'
    }
    
    # Pool the proportions across imputations using Rubin's Rules
    pooled_mean <- mean(prop_values)
    between_var <- var(prop_values)
    within_var <- mean(prop_variances)
    pooled_var <- within_var + (1 + 1 / length(prop_values)) * between_var
    pooled_se <- sqrt(pooled_var)
    
    # Update the hbdata template with the pooled proportion
    hbdata_new$value[hbdata_new$body_part == body_parts[i]] <- pooled_mean
    hbdata_new$se[hbdata_new$body_part == body_parts[i]] <- pooled_se
  }
  
  # Save the updated dataframe for this depression level
  hbdata_list[[depr_lvl]] <- hbdata_new
}

# Combine all dataframes for plotting
human_body_data <- bind_rows(
  hbdata_list %>%
    purrr::imap(~ .x %>% mutate(DEPR_LVL = .y))
)

# Plot faceted heatmaps with proportions
ggplot(human_body_data) +
  geom_sf(aes(fill = value), color = "black") +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  ) +
  labs(fill = "Proportion") +
  facet_wrap(~DEPR_LVL)



## TO plot line graph of proportion by DEPR_LVL 

# Template data
body_parts <- c("Right Arm", "Right Chest", "Neck", "Upper Sternum", 
                "Lower Sternum", "Left Chest", "Left Arm", "Epigastric")
variables <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", 
               "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")


# Initialize a data frame to store results for plotting
plot_data <- data.frame(
  DEPR_LVL = integer(),
  body_part = character(),
  proportion = numeric(),
  se = numeric()
)

# Iterate over depression levels (1 to 5)
for (depr_lvl in 1:5) {
  # Iterate over body region variables (CDQ009A to CDQ009H)
  for (i in seq_along(variables)) {
    # Placeholder for proportions and variances across imputations
    prop_values <- numeric()
    prop_variances <- numeric()
    
    # Process each imputed dataset
    dat_mice <- complete(imputations, "all")
    for (dat in dat_mice) {
      # Define the survey design for this dataset
      nhanes_design <- svydesign(
        id = ~SDMVPSU,
        weights = ~MEC15YR,
        strata = ~SDMVSTRA,
        data = dat,
        nest = TRUE
      )
      #Subset for those who responded to questionnaire of CP
      svy_subset1 <- subset(nhanes_design, !apply(
        nhanes_design$variables[, responses], 1,
        function(x) all(x == 0 | is.na(x))
      ))
      
      # Subset the survey data to the current depression level
      svy_subset <- subset(svy_subset1, DEPR_LVL == depr_lvl)
      
      # Calculate survey-weighted proportion for the current body region variable (e.g., CDQ009A, CDQ009B)
      prop <- svymean(as.formula(paste0("~", variables[i])), svy_subset, na.rm = TRUE)
      
      # Extract the proportion for '1' (second row)
      prop_values <- c(prop_values, as.numeric(prop[2]))  # prop[2] gives the proportion for '1'
      prop_variances <- c(prop_variances, attr(prop, "var")[2])  # Standard error for '1'
    }
    
    # Pool the proportions across imputations using Rubin's Rules
    pooled_mean <- mean(prop_values)
    between_var <- var(prop_values)
    within_var <- mean(prop_variances)
    pooled_var <- within_var + (1 + 1 / length(prop_values)) * between_var
    pooled_se <- sqrt(pooled_var)
    
    # Store the results in the plot_data dataframe
    plot_data <- rbind(plot_data, data.frame(
      DEPR_LVL = depr_lvl,
      body_part = body_parts[i],
      proportion = pooled_mean,
      se = pooled_se
    ))
  }
}

# Plot the line graph
ggplot(plot_data, aes(x = DEPR_LVL, y = proportion, color = body_part, group = body_part)) +
  geom_line(size = 1.2) +
  geom_point() +
  scale_color_manual(values = RColorBrewer::brewer.pal(8, "Set2")) +  # Color palette for 8 lines
  theme_minimal() +
  labs(
    title = "Proportions of Chest Pain Locations Across Depression Levels",
    x = "Depression Level",
    y = "Proportion",
    color = "Body Part"
  ) +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16)
  )



## Optional smooth plot

# Plot a smooth relationship between DEPR_LVL and proportion (separate plot)
smooth_plot <- ggplot(plot_data, aes(x = DEPR_LVL, y = proportion, color = body_part)) +
  geom_smooth(method = "loess", size = 1.2, se = FALSE) +  # Smoothed line, no confidence interval shading
  scale_color_manual(values = RColorBrewer::brewer.pal(8, "Set2")) +  # Color palette for 8 lines
  theme_minimal() +
  labs(
    title = "Smoothed Relationship of Proportion by Depression Level",
    x = "Depression Level",
    y = "Proportion",
    color = "Body Part"
  ) +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16)
  )