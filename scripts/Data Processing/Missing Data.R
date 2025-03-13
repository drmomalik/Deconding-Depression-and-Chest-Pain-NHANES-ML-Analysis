# Load necessary libraries
library(dplyr)        # Data manipulation
library(ggplot2)      # Data visualization
library(DataExplorer) # EDA utilities
library(psych)        # For descriptive statistics
library(here)

data <- subset_new

# Load your dataset (replace 'your_data.csv' with your file path)
data <- read.csv(here("your_data.csv"), header = TRUE)

# 1. View the first few rows of the dataset
print("First few rows of the dataset:")
head(data)

# 2. Check the structure of the dataset (data types, column names, etc.)
print("Structure of the dataset:")
str(data)

# 3. Summary statistics for each variable
print("Summary statistics:")
summary(data)

# 4. Descriptive statistics (mean, median, standard deviation, etc.)
print("Descriptive statistics:")
describe(data)

# 5. Check for missing values
print("Count of missing values per column:")
colSums(is.na(data))

# 6. Visualize the missing data pattern (if any) with percentages
plot_missing <- function(data) {
  missing_data <- colSums(is.na(data))
  total_rows <- nrow(data)
  missing_data <- missing_data[missing_data > 0]
  
  if (length(missing_data) > 0) {
    percentages <- round((missing_data / total_rows) * 100, 1)
    
    # Adjust margins to accommodate angled x-axis labels
    par(mar = c(12, 4, 4, 2))  # Bottom, Left, Top, Right margins
    
    # Set ylim to ensure space above the tallest bar
    ylim_max <- max(missing_data) * 1.2
    bar_positions <- barplot(
      missing_data, 
      main = "Missing Data Pattern", 
      ylab = "Count of Missing Values", 
      xlab = "",  # Remove x-axis label
      col = "skyblue", 
      ylim = c(0, ylim_max), # Extend y-axis range
      names.arg = rep("", length(missing_data)) # Hide default labels
    )
    
    # Add angled x-axis labels
    text(
      x = bar_positions, 
      y = par("usr")[3] - 0.02 * diff(par("usr")[3:4]), # Position slightly below the axis
      labels = names(missing_data), 
      srt = 45,  # Rotate text by 45 degrees
      adj = 1,   # Right-align text
      xpd = TRUE, # Allow drawing outside plot region
      cex = 0.8  # Smaller text size
    )
    
    # Overlay text showing the percentage of missing data
    text(
      x = bar_positions, 
      y = missing_data, 
      labels = paste0(percentages, "%"), 
      pos = 3, 
      cex = 0.8, 
      col = "red"
    )
  } else {
    print("No missing values found.")
  }
}




# Example usage:
plot_missing(data)


# 7. Visualize the distribution of numeric variables
numeric_cols <- select_if(data, is.numeric)
if (ncol(numeric_cols) > 0) {
  print("Histograms of numeric variables:")
  for (col in colnames(numeric_cols)) {
    print(ggplot(data, aes_string(col)) +
            geom_histogram(bins = 30, fill = "skyblue", color = "black") +
            ggtitle(paste("Distribution of", col)))
  }
} else {
  print("No numeric columns to plot.")
}

# 8. Box plots to detect outliers in numeric variables
print("Box plots for numeric variables to detect outliers:")
for (col in colnames(numeric_cols)) {
  print(ggplot(data, aes_string(x = "1", y = col)) +
          geom_boxplot(fill = "lightgreen", color = "black") +
          ggtitle(paste("Box plot of", col)) +
          xlab("") +
          theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()))
}

# 9. Distribution of categorical variables
categorical_cols <- select_if(data, is.factor)
if (ncol(categorical_cols) > 0) {
  print("Bar plots of categorical variables:")
  for (col in colnames(categorical_cols)) {
    print(ggplot(data, aes_string(col)) +
            geom_bar(fill = "orange", color = "black") +
            ggtitle(paste("Distribution of", col)))
  }
} else {
  print("No categorical columns to plot.")
}

# 10. Correlation matrix for numeric variables
if (ncol(numeric_cols) > 1) {
  print("Correlation matrix of numeric variables:")
  cor_matrix <- cor(numeric_cols, use = "pairwise.complete.obs")
  print(cor_matrix)
  # Optional: Visualize the correlation matrix
  library(reshape2)
  library(ggcorrplot)
  ggcorrplot(cor_matrix, hc.order = TRUE, type = "lower",
             lab = TRUE, lab_size = 3, method = "square",
             colors = c("blue", "white", "red"))
} else {
  print("Not enough numeric columns to compute a correlation matrix.")
}

# 11. Pairwise scatter plots for a quick view of relationships between variables
if (ncol(numeric_cols) > 1) {
  print("Pairwise scatter plots for numeric variables:")
  pairs(numeric_cols, main = "Pairwise Scatter Plots of Numeric Variables")
}

# Additional EDA Tools (Optional)
# DataExplorer::create_report(data) # Automatically generates an EDA report in HTML
