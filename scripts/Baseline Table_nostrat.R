### Create a table 1 without any stratification 


create_summary_table <- function(data, variables, var_names = NULL, level_names = NULL) {
  # Recode factor levels based on the level_names argument
  if (!is.null(level_names)) {
    for (var in names(level_names)) {
      if (var %in% colnames(data)) {
        data[[var]] <- factor(data[[var]], levels = names(level_names[[var]]), labels = level_names[[var]])
      }
    }
  }
  
  # Initialize a list to store rows for the table
  table_rows <- list()
  
  # Loop through variables to calculate values for the table
  for (var in variables) {
    var_name <- ifelse(!is.null(var_names) && var %in% names(var_names), var_names[[var]], var)
    
    # Handle categorical variables
    if (is.factor(data[[var]])) {
      counts <- table(data[[var]])
      proportions <- round(prop.table(counts) * 100, 1)
      
      # Add a row for the categorical variable name
      table_rows[[length(table_rows) + 1]] <- data.frame(
        Variable = var_name,
        Level = NA,
        Summary = NA
      )
      
      # Add rows for each level of the categorical variable
      for (level_name in names(counts)) {
        table_rows[[length(table_rows) + 1]] <- data.frame(
          Variable = NA,
          Level = level_name,
          Summary = paste0(counts[[level_name]], " (", proportions[[level_name]], "%)")
        )
      }
    }
    
    # Handle continuous variables
    else if (is.numeric(data[[var]])) {
      # Check for normality
      normal_test <- shapiro.test(data[[var]])$p.value
      if (normal_test > 0.05) {
        # Normally distributed: Use mean and SD
        mean_val <- round(mean(data[[var]], na.rm = TRUE), 1)
        sd_val <- round(sd(data[[var]], na.rm = TRUE), 1)
        summary_stat <- paste0(mean_val, " (SD ", sd_val, ")")
      } else {
        # Non-normally distributed: Use median and IQR
        median_val <- round(median(data[[var]], na.rm = TRUE), 1)
        iqr_val <- round(IQR(data[[var]], na.rm = TRUE), 1)
        summary_stat <- paste0(median_val, " [IQR: ", iqr_val, "]")
      }
      
      # Add row for continuous variable
      table_rows[[length(table_rows) + 1]] <- data.frame(
        Variable = var_name,
        Level = NA,
        Summary = summary_stat
      )
    }
  }
  
  # Combine rows into a data frame
  summary_table <- bind_rows(table_rows) %>%
    select(Variable, Level, Summary)
  
  return(summary_table)
}

# Example usage
level_names <- list(
  "urgency" = c("0" = "Elective", "1" = "Urgent"),
  "approach" = c("0" = "Sternotomy", "1" = "Mini-thoracotomy", "2" = "Robotic"),
  "annuloplasty" = c("0" = "Ring", "1" = "Band"),
  "technique_leaflet" = c("0" = "Neochords", "1" = "Resection", "2" = "Mixed"),
  "intraop_residual_mr" = c("0" = "None", "1" = "Mild", "2" = "Moderate")
)

summary_table <- create_summary_table(
  data = data_filt,
  variables = c("urgency", "approach", "ring_size", "total_pump", "xc_time", "intraop_residual_mr"),
  var_names = c(
    "urgency" = "Urgency",
    "approach" = "Surgical Approach",
    "ring_size" = "Ring Size",
    "total_pump" = "Total Pump Time",
    "xc_time" = "Cross Clamp Time",
    "intraop_residual_mr" = "Post-repair Residual MR"
  ),
  level_names = level_names
)

# Formatting the table using flextable
formatted_summary_table <- summary_table %>%
  flextable() %>%
  set_header_labels(
    Variable = "Covariate",
    Level = "",
    Summary = "Summary Statistics"
  ) %>%
  theme_vanilla() %>%
  autofit() %>%
  bold(part = "header") %>%
  bold(j = "Variable", part = "all") %>%  
  italic(j = "Level", part = "all") %>%
  align(j = c("Variable", "Level"), align = "left", part = "all") %>%
  align(j = "Summary", align = "center", part = "all")

# Displaying the formatted table
formatted_summary_table

# Export as image or Word document
save_as_docx(formatted_summary_table, path = here("summary_table.docx"))
save_as_image(formatted_summary_table, path = here("summary_table.png"))
