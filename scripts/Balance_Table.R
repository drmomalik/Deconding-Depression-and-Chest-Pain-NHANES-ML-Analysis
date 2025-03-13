library(flextable)
library(dplyr)


df <- final_balance
# Split Variable into "Variable Name" and "Level" (if applicable)
df <- df %>%
  separate(Variable, into = c("Variable_Name", "Variable_Level"), sep = ":", fill = "right") %>%
  mutate(Variable_Level = ifelse(is.na(Variable_Level), "", Variable_Level))  # Replace NA with empty string

# Reshape data for flextable
df_long <- df %>%
  rename(SMD_Population = SMD_before, 
         SMD_Matched = SMD_after) %>%
  select(Variable_Name, Variable_Level, 
         Treatment_Mean_before, Control_Mean_before, SMD_Population, 
         Treatment_Mean_after, Control_Mean_after, SMD_Matched)

# Rename columns for final table
colnames(df_long) <- c("Variable", "Level", 
                       "Depressed_Population", "Non_Depressed_Population", "SMD_Population", 
                       "Depressed_Matched", "Non_Depressed_Matched", "SMD_Matched")

# Create flextable with merged headers
flextable(df_long) %>%
  set_header_labels(
    Variable = "Variable Name",
    Level = "Variable Level",
    Depressed_Population = "Depressed",
    Non_Depressed_Population = "Non-Depressed",
    SMD_Population = "SMD",
    Depressed_Matched = "Depressed",
    Non_Depressed_Matched = "Non-Depressed",
    SMD_Matched = "SMD"
  ) %>%
  add_header(
    Depressed_Population = "Population Cohorts",
    Non_Depressed_Population = "Population Cohorts",
    SMD_Population = "Population Cohorts",
    Depressed_Matched = "Matched Cohorts",
    Non_Depressed_Matched = "Matched Cohorts",
    SMD_Matched = "Matched Cohorts",
    top = TRUE
  ) %>%
  merge_h(part = "header") %>%  # Merge headers for clarity
  theme_booktabs() %>%  # Apply clean table theme
  flextable::autofit()  # Adjust column widths automatically



flextable::