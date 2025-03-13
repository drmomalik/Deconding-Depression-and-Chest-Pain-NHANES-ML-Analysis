### Create a table 1 without any stratification 


create_summary_table <- function(svy_design, variables, var_names = NULL, level_names = NULL) {
  # Recode factor levels based on the level_names argument
  if (!is.null(level_names)) {
    for (var in names(level_names)) {
      if (var %in% colnames(svy_design$variables)) {
        svy_design$variables[[var]] <- factor(svy_design$variables[[var]], levels = names(level_names[[var]]), labels = level_names[[var]])
      }
    }
  }
  
  # Initialize a list to store rows for the table
  table_rows <- list()
  
  # Loop through variables to calculate values for the table
  for (var in variables) {
    var_name <- ifelse(!is.null(var_names) && var %in% names(var_names), var_names[[var]], var)
    
    # Handle categorical variables
    if (is.factor(svy_design$variables[[var]])) {
      counts <- table(svy_design$variables[[var]])
      proportions <- round(prop.table(counts) * 100, 1)
      
      prop_table <- svytable(as.formula(paste("~", var)), svy_design)
      proportions1 <- round(prop_table / sum(prop_table) * 100, 1)
      levels <- names(prop_table)
      
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
          Summary = paste0(counts[[level_name]], " (", proportions[[level_name]], "%)"),
          svy.wt = paste0(proportions1[level_name], "%")
        )
      }
    }
    
    # Handle continuous variables
    else if (is.numeric(svy_design$variables[[var]])) {
        # Normally distributed: Use mean and SD
        mean_val <- round(mean(svy_design$variables[[var]], na.rm = TRUE), 1)
        sd_val <- round(sd(svy_design$variables[[var]], na.rm = TRUE), 1)
        summary_stat <- paste0(mean_val, " (SD ", sd_val, ")")
      
      mean_val <- svymean(as.formula(paste("~", var)), svy_design, na.rm = TRUE)
      mean_stat <- round(coef(mean_val), 1)
      se_stat <- round(SE(mean_val), 1)
      
      # Add row for continuous variable
      table_rows[[length(table_rows) + 1]] <- data.frame(
        Variable = var_name,
        Level = NA,
        Summary = summary_stat,
        svy.wt = paste0(mean_stat, " (", se_stat, ")")
      )
    }
  }
  
  # Combine rows into a data frame
  summary_table <- bind_rows(table_rows) %>%
    select(Variable, Level, Summary, svy.wt)
  
  return(summary_table)
}

# Example usage
level_names <- list(
  "DMDBORNT" = c("1" = "USA", "2" = "Other"),
  "RIAGENDR" = c("1" = "Female", "2" = "Male"),
  "DMDEDUC2" = c("1" = "< 9th grade", "2" = "9-11th grade (no diploma)", "3" = "High School Grad/GED", "4" = "Some College or AA", "5" = "College Graduate"),
  "DEPR_LVL" = c("1" = "None/Mild (0-4)", "2" = "Mild (5-9)", "3" = "Moderate (10-14)", "4" = "Moderately Severe (15-19)", "5" = "Severe (20-27)"),
  "AGE_BIN" = c("1" = "40-49 years", "2" = "50-59 years", "3" = "60-69 years", "4" = "70-79 years", "5" = "80 years or older"),
  "PAQMV" = c("0" = "Less than weekly", "1" = "At least weekly"),
  "RIDRETH1" = c("1" = "Mexican American", "2" = "Other Hispanic", "3" = "Non-Hispanic White", "4" = "Non-Hispanic Black", "5" = "Other Race"),
  "INC3" = c("1" = "0-0.99", "2" = "1-1.99", "3" = "2 or greater"),
  "HIQ011" = c("1" = "Yes", "2" = "No"),
  "DEPR_BIN" = c("1" = "PHQ-9 Score < 15", "2" = "PHQ-9 Score  ≥ 15"),
  "SDDSRVYR" = c("4" = "2005-2006", "5" = "2007-2008", "6" = "2009-2010", "7" = "2011-2012", "8" = "2013-2014", "9" = "2015-2016", "66" = "2017-2020 (Pre-pandemic)"),
  "MEDDEP" = c("0" = "No", "1" = "Yes"),
  "DIDTOT" = c("0" = "No", "1" = "Borderline", "2" = "Diet-Controlled", "3" = "Medication-Controlled", "4" = "Insulin-dependent"),
  "BPQ020" = c("1" = "Yes", "2" = "No"),
  "BPQ080" = c("1" = "Yes", "2" = "No"),
  "SMQ020" = c("1" = "Yes", "2" = "No"),
  "SMQ040" = c("1" = "Everyday", "2" = "Some days", "3" = "No"),
  "DUQTOT" = c("1" = "Yes", "0" = "No")
)

summary_table <- create_summary_table(
  svy_design = svy_design,
  variables = c("AGE_BIN", "RIAGENDR", "BMXBMI", "PAQMV", "RIDRETH1", "DMDBORNT" , "DMDEDUC2", "INC3", "HIQ011", 
                "DEPR_LVL", "DEPR_TOT" , "DEPR_BIN", "SDDSRVYR", "MEDDEP", "DIDTOT", "BPQ020", 
                "BPQ080", "SMQ020", "SMQ040", "ALQ130", "DUQTOT"),
  var_names = c("AGE_BIN" = "Age", "RIAGENDR" = "Sex", "PAQMV" = "Physical Activity (Moderate-Vigorous)",
                "RIDRETH1" = "Race/Ethnicity", "DMDBORNT" = "Country of Birth", "DMDEDUC2" = "Education Level", 
                "INC3" = "Income (as a proportion of the poverty threshold)", "HIQ011" = "Health Insurance Coverage", "SDDSRVYR" = "NHANES Cycle", 
                "DEPR_TOT" = "PHQ-9 Score (Total)", "DEPR_BIN" = "PHQ-9 Score (Binary)", "MEDDEP" = "Antidepressant-Use",
                "DIDTOT" = "Diabetes", "BPQ020" = "Hypertension", "BPQ080" = "Hypercholesterolemia", "SMQ020" = "Smoked at least 100 Cigarettes", 
                "SMQ040" = "Currently Smoking", "ALQ130" = "Avg. # of alcholic drinks per day over last year", "DUQTOT" = "Any history of illicit drug use",
                "DEPR_LVL" = "Depression Severity (based on PHQ-9 Score)",
                "BMXBMI" = "BMI (kg/m²)"
  ),
  level_names = level_names
)


design <- svydesign(
  id = ~SDMVPSU,
  weights = ~MEC15YR,
  strata = ~SDMVSTRA,
  data = wd,
  nest = TRUE
)

svy_design <- subset(design, !apply(
  design$variables[, responses], 1,
  function(x) all(x == 0 | is.na(x))
))


# Formatting the table using flextable
formatted_summary_table <- summary_table %>%
  flextable() %>%
  set_header_labels(
    Variable = "Covariate",
    Level = "",
    Summary = "Total Cohort (n = 1880)",
    svy.wt = "Survey Weighted\nProportions/Means"
  ) %>%
  theme_vanilla() %>%
  autofit() %>%
  bold(part = "header") %>%
  bold(j = "Variable", part = "all") %>%  
  italic(j = "Level", part = "all") %>%
  align(j = c("Variable", "Level"), align = "left", part = "all") %>%
  align(j = "Summary", align = "center", part = "all") %>% 
  align(j = "svy.wt", align = "center", part = "all")


# Displaying the formatted table
formatted_summary_table

# Export as image or Word document
save_as_docx(formatted_summary_table, path = here("summary_table.docx"))
save_as_image(formatted_summary_table, path = here("image/summary_table_svy.png"))





##### 


## Create svydesign object
design <- svydesign(
  id = ~SDMVPSU,
  weights = ~MEC15YR,
  strata = ~SDMVSTRA,
  data = wd,
  nest = TRUE
)

svy_design <- subset(design, !apply(
  design$variables[, responses], 1,
  function(x) all(x == 0 | is.na(x))
))


# Change tot_reg to continuous variable
svy_design$variables$TOT_REG <- as.factor(svy_design$variables$TOT_REG)


# Example usage
level_names <- list(
  "CDQ009A" = c("1" = "Yes", "0" = "No"),
  "CDQ009B" = c("1" = "Yes", "0" = "No"),
  "CDQ009C" = c("1" = "Yes", "0" = "No"),
  "CDQ009D" = c("1" = "Yes", "0" = "No"),
  "CDQ009E" = c("1" = "Yes", "0" = "No"),
  "CDQ009F" = c("1" = "Yes", "0" = "No"),
  "CDQ009G" = c("1" = "Yes", "0" = "No"),
  "CDQ009H" = c("1" = "Yes", "0" = "No"),
  "CDQ008" = c("1" = "Yes", "0" = "No"),
  "CDQ010" = c("1" = "Yes", "0" = "No"),
  "ROSE" = c("0" = "No", "1" = "Grade 1", "2" = "Grade 2"),
  "CAD" = c("0" = "No", "1" = "Yes"),
  "CADTOT" = c("0" = "No", "1" = "Yes")
)


summary_table2 <- create_summary_table(
  svy_design = svy_design,
  variables = c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H", "TOT_REG", 
                "CDQ008" , "CDQ010", "ROSE", "CAD", "CADTOT"),
  var_names = c("CDQ009A" = "Right Arm", "CDQ009B" = "Right Chest", "CDQ009C" = "Neck", "CDQ009D" = "Upper Sternum", 
                "CDQ009E" = "Lower Sternum", "CDQ009F" = "Left Chest", "CDQ009G" = "Left Arm", "CDQ009H" = "Epigastric", 
                "TOT_REG" = "Total Regions of Chest Pain", 
                "CDQ008" = "Severe chest pain lasting more than 30 minutes", "CDQ010" = "Shortness of breath on stairs/incline",
                "ROSE" = "Rose Criteria Angina", "CAD" = "History of CAD/MI/Angina", 
                "CADTOT" = "History of CAD/MI/Angina and/or Rose Criteria Grade 1-2 Angina"
  ),
  level_names = level_names
)

# Formatting the table using flextable
formatted_summary_table <- summary_table2 %>%
  flextable() %>%
  set_header_labels(
    Variable = "Chest Pain Questionnaire responses",
    Level = "",
    Summary = "Total Cohort (n = 1880)",
    svy.wt = "Survey Weighted\nProportions"
  ) %>%
  theme_vanilla() %>%
  autofit() %>%
  bold(part = "header") %>%
  bold(j = "Variable", part = "all") %>%  
  italic(j = "Level", part = "all") %>%
  align(j = c("Variable", "Level"), align = "left", part = "all") %>%
  align(j = "Summary", align = "center", part = "all") %>% 
  align(j = "svy.wt", align = "center", part = "all")


# Displaying the formatted table
formatted_summary_table

save_as_image(formatted_summary_table, path = here("image/summary_table2_svy.png"))
