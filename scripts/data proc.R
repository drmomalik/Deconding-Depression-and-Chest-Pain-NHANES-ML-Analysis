library(tidyverse)
library(haven)
library(here)

dpq_J <- nhanes("DPQ_J", translated = F) # for year 2017-2018

dpq_vars <- c("DPQ010", "DPQ020", "DPQ030", "DPQ040", "DPQ050", 
              "DPQ060", "DPQ070", "DPQ080", "DPQ090")
# Filter out rows where any of the DPQ variables are 7 or 9 using package 'dplyr' (refused/don't know)
df_fin <- dpq_J %>%
  dplyr::mutate(across(all_of(dpq_vars), ~na_if(., 7))) %>%
  dplyr::mutate(across(all_of(dpq_vars), ~na_if(., 9)))

###### variable transformation

df_fin <- df_fin %>% mutate(DEPR_TOT = DPQ010 + DPQ020 + DPQ030 + DPQ040 + DPQ050 + DPQ060 + DPQ070 + DPQ080 + DPQ090)

### Create PHQ-9 bins
df_fin$DEPR_LVL <- cut(df_fin$DEPR_TOT,
                       breaks = c(-Inf, 4, 9, 14, 19, Inf),
                       labels = 1:5,
                       right = TRUE) 


## Import outcome variable of chest pain presence and location
df_cp <- nhanes("CDQ_J", translated = F)
df_cp <- c("CDQ001", "CDQ002", "CDQ003",  "CDQ004",  "CDQ005",  "CDQ006",  "CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E",
           "CDQ009F", "CDQ009G", "CDQ009H", "CDQ008", "CDQ010")
df_cp <- df_cp %>%
  dplyr::mutate(across(all_of(cp_var), ~na_if(., 77))) %>%
  dplyr::mutate(across(all_of(cp_var), ~na_if(., 99)))

## Merge data frames
merge_df <- merge(df_cp, wd_na, by = "SEQN", all = FALSE)

##Keep only those that have had chest pain
df <- subset(merge_df, CDQ001 == 1)

# Subset dataframe with participants that responded to any of the pain location questions
df_lim <- df[rowSums(is.na(df[, c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", 
                                  "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")])) < 8, ]


### Repeat above steps but for all years between 2005 - 2020 (Prepandemic)

# Define variables and years
vars <- c("DPQ_", "CDQ_", "DEMO_", "DIQ_",
          "PAQ_", "SMQ_", "BPQ_", "HIQ_",
          "BMX_", "INQ_", "MCQ_", "HUQ_",
          "ALQ_", "DUQ_", "RXQ_RX_"
)
years <- c("D", "E", "F", "G", "H", "I")

codes <- c("a70010", "c00249", "c00250", "d00144", "d00145", "d00146", "d00181", "d00217", 
           "d00236", "d00259", "d00395", "d00873", "d00874", "d00875", "d00876", "d00877", 
           "d00880", "d00882", "d00883", "d00884", "d00976", "d03157", "d03181", "d03804", 
           "d03808", "d04025", "d04332", "d04812", "d05355", "d06635", "d07113", "d07740", 
           "d08114", "d08125", "h00035")

# Initialize an empty list to hold yearly combined datasets
yearly_combined_list <- list()

# Loop through years
for (year in years) {
  year_data_list <- list() # Temp storage for variable data for this year
  
  # Loop through variables
  for (var in vars) {
    var_year <- paste0(var, year) # Combine variable and year
    
    message("Fetching data for: ", var_year)
    
    # Fetch data
    data <- tryCatch(
      read_xpt(here(paste0("data/",var,year,".xpt"))),
      error = function(e) {
        message("Error fetching ", var_year, ": ", e$message)
        return(NULL) # Skip if error occurs
      }
    )
    # Extract the use of anti-depressant 
    if (var == "RXQ_RX_") {
      # Create MEDDEP column
      data$MEDDEP <- ifelse(data$RXDDRGID %in% codes, 1, 0)
      
      # select only SEQN and MEDDEP 
      data <- data %>% dplyr::select(SEQN, MEDDEP)
      
      #collapse over unqiue SEQN 
      data <- data %>%
        dplyr::group_by(SEQN) %>%               # Group by SEQN
        dplyr::summarize(MEDDEP = max(MEDDEP), .groups = "drop") # Take the maximum value of MEDDEP for each SEQN
    }
    
    # Append to list if data exists
    if (!is.null(data)) {
      year_data_list[[var]] <- data
    } else {
      print("Data frame is empty!")
    }
  }
  
  # Merge all variable datasets for this year by SEQN
  if (length(year_data_list) > 0) {
    combined_year_data <- reduce(year_data_list, full_join, by = "SEQN")
    combined_year_data <- mutate(combined_year_data, YEAR = year) # Add year identifier
    yearly_combined_list[[year]] <- combined_year_data
  }
}

# Merge across cycles 
# Combine all yearly datasets into one, aligning by columns
final_data <- bind_rows(yearly_combined_list)



## Add prepandemic 2019-2020 data
# Define the new year prefix and variables
new_year_prefix <- "P_"
vars <- c("DPQ", "CDQ", "DEMO", "DIQ",
          "PAQ", "SMQ", "BPQ", "HIQ", "HUQ",
          "BMX", "CRP", "INQ", "MCQ", 
          "ALQ", "DUQ", "RXQ_RX")

# Initialize a list to store new data
new_year_data_list <- list()

# Loop through variables for the new year
for (var in vars) {
  var_name <- paste0(new_year_prefix, var) # Construct variable name
  message("Fetching data for: ", var_name)
  
  # Fetch data
  data <- tryCatch(
    read_xpt(here(paste0("data/",var_name,".xpt"))),
    error = function(e) {
      message("Error fetching ", var, ": ", e$message)
      return(NULL) # Skip if error occurs
    }
  )
  # Extract the use of anti-depressant 
  if (var == "RXQ_RX") {
    # Create MEDDEP column
    data$MEDDEP <- ifelse(data$RXDDRGID %in% codes, 1, 0)
    
    # select only SEQN and MEDDEP 
    data <- data %>% dplyr::select(SEQN, MEDDEP)
    
    #collapse over unqiue SEQN 
    data <- data %>%
      dplyr::group_by(SEQN) %>%               # Group by SEQN
      dplyr::summarize(MEDDEP = max(MEDDEP), .groups = "drop") # Take the maximum value of MEDDEP for each SEQN
  }
  
  # Append to list if data exists
  if (!is.null(data)) {
    new_year_data_list[[var]] <- data
  } else {
    print("Data is empty!")
  }
}

# Combine all variables for the new year into a single dataframe
new_year_data <- reduce(new_year_data_list, full_join, by = "SEQN")
new_year_data <- mutate(new_year_data, YEAR = "P") # Add year identifier

# Append new year data to the existing dataframe
final_data <- bind_rows(final_data, new_year_data)

# Verify the structure of the updated dataframe
str(final_data)

## selecting relevant variables for data frame
# List of variables to include in the new dataframe
selected_vars <- c(
  "RIAGENDR", "RIDAGEYR", "RIDRETH1", "DMDBORN", "DMDBORN2",
  "DMDBORN4", "DMDEDUC2", "DMDMARTL",
  "WTMECPRP", "WTMEC2YR", "SDMVPSU", "INDFMPIR", "DIQ010", "DID040", 
  "DIQ050", "DID060", "DIQ070", "PAQ650", "PAQ665", "SMQ020", "SMD030", "SMQ040", "SMD650", 
  "IND235", "INDFMMPC", "HIQ011", "BPQ020", "BPQ040A", "BMXBMI", 
  "BPQ080", "MCQ160C", "MCQ180C", "MCQ160D", "MCQ180D", "MCQ160E", "MCQ180E", 
  "DPQ010", "DPQ020", "DPQ030", "DPQ040", "DPQ050", "DPQ060", "DPQ070", 
  "DPQ080", "DPQ090", "CDQ001", "CDQ002", "CDQ003", "CDQ004", "CDQ005", 
  "CDQ006", "CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", 
  "CDQ009F", "CDQ009G", "CDQ009H", "CDQ008", "CDQ010", "SEQN", "SDDSRVYR",
  "ALQ130", "DUQ200", "DUQ250", "DUQ290", "DUQ330", "DUQ370", "MEDDEP", "SDMVSTRA",
  "PAD200", "PAD320", "HUQ090"
  )

# Create a new dataframe with only the selected variables
df_lim <- final_data[, selected_vars, drop = FALSE]

### If needed, to check if selected_vars are columns in the final data
missing_vars <- setdiff(selected_vars, colnames(final_data))









## set new df name for filtered "dont know" and "missing"
wd_na <- df_lim

### Optional for removing rows with 7,9 or 77,99 in DPQ and CDQ
miss1 <- c("DIQ010", "DPQ010", "DPQ020", "DPQ030", "DPQ040", "DPQ050", 
           "DPQ060", "DPQ070", "DPQ080", "DPQ090", "INDFMMPC",
           "SMQ020", "HIQ011", "BPQ080", "DUQ200", "DUQ250",
           "DUQ370", "DMDBORN2","DMDBORN", "SMQ040", "BPQ020",
           "CDQ001", "CDQ002", "CDQ003",  "CDQ004",  "CDQ005",  "CDQ006",
           "DUQ290" , "DUQ330", "DMDEDUC2", "PAD200", "PAD320", "HUQ090")
miss2 <- c("CDQ001", "CDQ002", "CDQ003",  "CDQ004",  "CDQ005",  "CDQ006",  "CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E",
           "CDQ009F", "CDQ009G", "CDQ009H", "CDQ008", "CDQ010", "DMDBORN4", "DMDMARTL")

# Filter out rows where any of the DPQ variables are 7 or 9 using package 'dplyr' (refused/don't know)
# Create a function to count changes to NA
count_na_changes <- function(data, vars, value_to_replace) {
  sum(data[vars] == value_to_replace, na.rm = TRUE)
}

# Count changes for the first block of replacements
miss1_na_7 <- count_na_changes(wd_na, miss1, 7)
miss2_na_77 <- count_na_changes(wd_na, miss2, 77)

# Apply the first set of replacements
wd_na <- wd_na %>%
  dplyr::mutate(across(all_of(miss1), ~na_if(., 7))) %>%
  dplyr::mutate(across(all_of(miss2), ~na_if(., 77)))

# Count changes for the second block of replacements
miss1_na_9 <- count_na_changes(wd_na, miss1, 9)
miss2_na_99 <- count_na_changes(wd_na, miss2, 99)

# Apply the second set of replacements
wd_na <- wd_na %>%
  dplyr::mutate(across(all_of(miss1), ~na_if(., 9))) %>%
  dplyr::mutate(across(all_of(miss2), ~na_if(., 99)))

# Total NA changes
total_na_changes <- miss1_na_7 + miss2_na_77 + miss1_na_9 + miss2_na_99

# Print the results
cat("Total cells changed to NA:\n",
    "- Miss 1 vars (7):", miss1_na_7, "\n",
    "- Miss 1 vars (77):", miss2_na_77, "\n",
    "- Miss 2 vars (9):", miss1_na_9, "\n",
    "- Miss 2 vars (99):", miss2_na_99, "\n",
    "Total:", total_na_changes, "\n")

###### variable transformation

## Create BMI bins 
wd_na$BMI_LVL <- cut(wd_na$BMXBMI,
                      breaks = c(-Inf, 18.4, 24.9, 29.9, 34.9, 39.9, Inf),
                      labels = 1:6,
                      right = TRUE)

wd_na <- wd_na %>% mutate(DEPR_TOT = DPQ010 + DPQ020 + DPQ030 + DPQ040 + DPQ050 + DPQ060 + DPQ070 + DPQ080 + DPQ090)

### Create PHQ-9 bins
wd_na$DEPR_LVL <- cut(wd_na$DEPR_TOT,
                      breaks = c(-Inf, 4, 9, 14, 19, Inf),
                      labels = 1:5,
                      right = TRUE) 
wd_na$DEPR_BIN <- cut(wd_na$DEPR_TOT,
                      breaks = c(-Inf, 14, Inf),
                      labels = 1:2,
                      right = TRUE) 



#Turn NA in chest pain location into 0 and presence of pain in location to 1
wd_na <- wd_na %>%
  mutate(across(c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", 
                  "CDQ009F", "CDQ009G", "CDQ009H"), 
                ~ ifelse(is.na(.) | . == 0, 0, 1)))

#Turn NA in chest pain location into 0 and presence of pain in location to 1
wd_na <- wd_na %>%
  mutate(across(c("CDQ008", "CDQ010"), 
                ~ ifelse(is.na(.) | . == 2, 0, 1)))


## Creating weighting variable for dataset (need to account for pre-pandemic dataset)
wd_na$MEC15YR <- ifelse(wd_na$SDDSRVYR == 66, wd_na$WTMECPRP * 3.2/15.2, wd_na$WTMEC2YR * 2/15.2)


# Combine DMDBORN, DMDBORN2, and DMDBORN4 into one column
wd_na$DMDBORN2[wd_na$DMDBORN2 %in% c(4, 5)] <- 2
wd_na$DMDBORN[wd_na$DMDBORN %in% 3] <- 2
wd_na$DMDBORNT <- coalesce(wd_na$DMDBORN, wd_na$DMDBORN2, wd_na$DMDBORN4)
wd_na <- wd_na %>% dplyr::select(-DMDBORN, -DMDBORN2, -DMDBORN4)

# Combine moderate and vigorous PA into one variable
wd_na$PAQMV <- ifelse(
  rowSums(sapply(c("PAQ650", "PAQ665", "PAD200", "PAD320"), function(x) wd_na[[x]] == 1), na.rm = TRUE) > 0,
  1,
  ifelse(
    rowSums(sapply(c("PAQ650", "PAQ665", "PAD200", "PAD320"), function(x) !is.na(wd_na[[x]])), na.rm = TRUE) == 0,
    NA,
    0
  )
)



# Create ROSE angina variable based on Rose Questionnaire criteria 
wd_na <- wd_na %>%
  mutate(
    ROSE = case_when(
      # Grade 1 Angina
      CDQ001 == 1 & CDQ002 == 1 & CDQ003 != 1 & CDQ004 == 1 & CDQ005 == 1 & CDQ006 == 1 &
        ((CDQ009D == 1 | CDQ009E == 1) | (CDQ009F == 1 & CDQ009G == 1)) ~ 1,
      
      # Grade 2 Angina
      CDQ001 == 1 & CDQ002 == 1 & CDQ003 == 1 & CDQ004 == 1 & CDQ005 == 1 & CDQ006 == 1 &
        ((CDQ009D == 1 | CDQ009E == 1) | (CDQ009F == 1 & CDQ009G == 1)) ~ 2,
      
      # Neither Grade 1 nor Grade 2 Angina
      TRUE ~ 0
    )
  )


# Create a new variable (CAD) for if respondent said yes to at least one of MCQ160C, MCQ160D, MCQ160E
wd_na$CAD <- ifelse(wd_na$MCQ160C == 1 | wd_na$MCQ160D == 1 | wd_na$MCQ160E == 1, 1, 0)

# Create variable if respondent positive for CAD or + ROSE questionnaire 
wd_na$CADTOT <- ifelse(wd_na$CAD == 1 | wd_na$ROSE == 2 | wd_na$ROSE == 1, 1, 0)


# Diabetic categorical variable - None, borderline, diet, medication, insulin +/- meds
wd_na$DIDTOT <- ifelse(wd_na$DIQ010 == 2,0, ifelse(wd_na$DIQ010 == 3, 1, ifelse(wd_na$DIQ050 == 1, 4, 
                                                                                ifelse(wd_na$DIQ070 == 1, 3, 2))))

# New definition for prevalence of depression
wd_na$DEPR_NEW <- ifelse(rowSums(is.na(wd_na[, c("DEPR_LVL", "MEDDEP", "HUQ090")])) == 3, NA, 
                         ifelse((wd_na$DEPR_TOT > 4 & wd_na$MEDDEP == 1) | (wd_na$DEPR_TOT > 4 & wd_na$HUQ090 == 1), 1, 0))

# Define atypical chest pain
wd_na$ATCP <- ifelse(rowSums(is.na(wd_na[, c("CDQ002", "CDQ003", "CDQ008", "CDQ009H", "CDQ009A", "CDQ009B")])) == 6, NA, 
                     ifelse(wd_na$CDQ002 == 2 | wd_na$CDQ003 == 2 | wd_na$CDQ008 == 1 | 
                              wd_na$CDQ009H == 1 | wd_na$CDQ009A == 1 | wd_na$CDQ009B == 1, 1, 0))

# Define typical chest pain
wd_na$TCP <- ifelse(rowSums(is.na(wd_na[, c("CDQ002", "CDQ003", "CDQ009D", "CDQ009E", "CDQ009F")])) == 5, NA, 
                    ifelse(wd_na$CDQ002 == 1 | (wd_na$CDQ003 == 1 & (wd_na$CDQ009D == 1 | 
                                                                       wd_na$CDQ009E == 1 | wd_na$CDQ009F == 1)), 1, 0))

  
# Creating drug use categorical variable (ever used any illicit drugs) 
wd_na$DUQTOT <- apply(wd_na[, c("DUQ200", "DUQ250", "DUQ290", "DUQ330", "DUQ370")], 1, function(row) {
  if (all(is.na(row))) {
    return(NA)  # All values are NA
  } else if (any(row == 1, na.rm = TRUE)) {
    return(1)  # At least one value is 1
  } else {
    return(0)  # All values are 0 or NA
  }
})


# Create bins for AGE (RIDAGEYR) (40-49, 50-59, 60-69, 70-79, >=80)
wd_na$AGE_BIN <- cut(wd_na$RIDAGEYR,
                      breaks = c(-Inf, 49, 59, 69, 79, Inf),
                      labels = 1:5,
                      right = TRUE) 

# Create bins for INDFMPIR (0-0.99, 1-1.99, 2-2.99, 3-3.99, 4-4.99, >5)
wd_na$INC_BIN <- cut(wd_na$INDFMPIR,
                     breaks = c(-Inf, 0.99, 1.99, 2.99, 3.99, 4.99, Inf),
                     labels = 1:6,
                     right = TRUE)

wd_na$INC3 <- cut(wd_na$INDFMPIR,
                  breaks = c(-Inf, 0.99, 1.99, Inf),
                  labels = 1:3,
                  right = TRUE)



# Create a variable to track total regions where pain is experienced
wd_na$TOT_REG <- rowSums(wd_na[,c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", 
                                  "CDQ009F", "CDQ009G", "CDQ009H")])


## create working dataset = wd
# remove unneeded variables 
wd <- wd_na %>% dplyr::select(-WTMECPRP, -WTMEC2YR, -PAQ650, -PAQ665, -BPQ040A, 
                        -MCQ160C, -MCQ180C, -MCQ160D, -MCQ180D, -MCQ160E, -MCQ180E,
                        -DPQ010, -DPQ020, -DPQ030, -DPQ040, -DPQ040, -DPQ050, -DPQ060,
                        -DPQ070, -DPQ080, -DPQ090, -DID040, -DID060, -SMD030
                        , -SMD650, -IND235, -DID040, -DID060, -DIQ050, -DIQ070, -SMD030
                        )


## setting correct class for each variable
# Function to reclassify variables
reclassify_variables <- function(data, factor_vars = NULL, numeric_vars = NULL) {
  for (var in factor_vars) {
    if (var %in% names(data)) {
      data[[var]] <- as.factor(data[[var]])
    } else {
      warning(paste("Variable", var, "not found in the dataset."))
    }
  }
  
  for (var in numeric_vars) {
    if (var %in% names(data)) {
      data[[var]] <- as.numeric(as.character(data[[var]])) # Ensure conversion from factors/characters
    } else {
      warning(paste("Variable", var, "not found in the dataset."))
    }
  }
  return(data)
}

# Define variables to reclassify
factor_vars <- c("RIAGENDR", "RIDRETH1",
                 "DIQ010", "SMQ020", "SMQ040",
                 "INDFMMPC", "HIQ011", "BPQ020",
                 "BPQ080", "CDQ001", "CDQ002",
                 "CDQ003", "CDQ004", "CDQ005",
                 "CDQ006", "CDQ009A", "CDQ009B",
                 "CDQ009C", "CDQ009D", "CDQ009E",
                 "CDQ009F", "CDQ009G", "CDQ009H", "DMDEDUC2", "DMDMARTL",
                 "CDQ008", "CDQ010", "SEQN", "SDDSRVYR",
                 "DUQ200", "DUQ250", "DUQ290", "DUQ330",
                 "DUQ370", "DUQTOT", "MEDDEP", "DEPR_LVL", "DMDBORNT",
                 "PAQMV", "CAD", "ROSE","DIDTOT", "AGE_BIN", "CADTOT", "ATCP",
                  "INC_BIN","SDMVPSU", "SDMVSTRA", "ALQ130",
                 "DEPR_NEW", "TCP", "HUQ090")  # Variables to convert to factors
numeric_vars <- c("RIDAGEYR", "INDFMPIR",
                  "BMXBMI", "ALQ130", "DEPR_TOT",
                  "MEC15YR", "TOT_REG")      # Variables to convert to numeric

# Reclassify variables in the dataframe
wd <- reclassify_variables(wd, factor_vars = factor_vars, numeric_vars = numeric_vars)



###  OPTIONAL: Recommended to subset the svydesign object

### Subsetting dataframe for only chest pain and those responding to location of CP
# Subsetting for those who responded to the chest pain questionnaire 
subsetd <- subset(wd, !is.na(CDQ001))

# Subsetting for those who have had chest pain
subsetd_1 <- subset(subsetd, CDQ001 == 1)

#Subset for chest pain and depression 
subset_new <- subset(subsetd, CDQ001 == 1 & DEPR_NEW == 1)
subset_new <- subset_new %>% select("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H",
                                    "CDQ008" , "CDQ010", "AGE_BIN", "RIAGENDR", "BMXBMI", "PAQMV", "RIDRETH1", "DMDBORNT" , 
                                    "DMDEDUC2", "INC3", "HIQ011", 
                                   "SDDSRVYR", "DIDTOT", "BPQ020", 
                                    "BPQ080", "SMQ020", "SMQ040", "ALQ130", "DUQTOT", "CAD")


# Subset dataframe with participants who responded to more than one pain location question
subsetd_2 <- subsetd_1[rowSums(subsetd_1[, c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", 
                                             "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")] != 0, na.rm = TRUE) > 0, ]



#Subset dataframe for MDA (only selecting for relevant variables)
subsetd_3 <- subsetd_2 %>% select("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H", "TOT_REG", 
                                 "CDQ008" , "CDQ010", "AGE_BIN", "RIAGENDR", "BMXBMI", "PAQMV", "RIDRETH1", "DMDBORNT" , "DMDEDUC2", "INC3", "HIQ011", 
                                 "DEPR_LVL", "DEPR_TOT" , "DEPR_BIN", "SDDSRVYR", "MEDDEP", "DIDTOT", "BPQ020", 
                                 "BPQ080", "SMQ020", "SMQ040", "ALQ130", "DUQTOT")

## Rename variables 

var_names = c("AGE_BIN" = "Age", "RIAGENDR" = "Sex", "PAQMV" = "Physical Activity (Moderate-Vigorous)",
              "RIDRETH1" = "Race/Ethnicity", "DMDBORNT" = "Country of Birth", "DMDEDUC2" = "Education Level", 
              "INC3" = "Income (as a proportion of the poverty threshold)", "HIQ011" = "Health Insurance Coverage", "SDDSRVYR" = "NHANES Cycle", 
              "DEPR_TOT" = "PHQ-9 Score (Total)", "DEPR_BIN" = "PHQ-9 Score (Binary)", "MEDDEP" = "Antidepressant-Use",
              "DIDTOT" = "Diabetes", "BPQ020" = "Hypertension", "BPQ080" = "Hypercholesterolemia", "SMQ020" = "Smoked at least 100 Cigarettes", 
              "SMQ040" = "Currently Smoking", "ALQ130" = "Avg. # of alcholic drinks per day over last year", "DUQTOT" = "Any history of illicit drug use",
              "DEPR_LVL" = "Depression Severity (based on PHQ-9 Score)",
              "BMXBMI" = "BMI (kg/m²)", "CDQ009A" = "Right Arm", "CDQ009B" = "Right Chest", "CDQ009C" = "Neck", "CDQ009D" = "Upper Sternum", 
              "CDQ009E" = "Lower Sternum", "CDQ009F" = "Left Chest", "CDQ009G" = "Left Arm", "CDQ009H" = "Epigastric", 
              "TOT_REG" = "Total Regions of Chest Pain", 
              "CDQ008" = "Severe chest pain lasting more than 30 minutes", "CDQ010" = "Shortness of breath on stairs/incline")
# Rename the columns in subsetd_3 using the var_names mapping
colnames(subsetd_3) <- ifelse(
  colnames(subsetd_3) %in% names(var_names), 
  var_names[colnames(subsetd_3)], 
  colnames(subsetd_3) # Keep columns not listed in var_names unchanged
)
###


