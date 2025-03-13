### OR Forest PLots for Simple and Full Model
library(tidyverse)
library(forestplot)

# Load necessary packages
library(ggplot2)
library(dplyr)

# Create a data frame for the forest plot with correct row names
model_data <- list(
  CDQ009A = data.frame(
    Variable = c("PHQ-9", "PHQ-9:AD", "PHQ-9:CAD", "Marginal Effect"),
    Mean = c(1.0227031, 1.0134347, 0.9983878, 1.0248963),
    Lower = c(0.9579024, 0.9379228, 0.9321012, 0.9609510),
    Upper = c(1.091888, 1.095026, 1.069388, 1.093097),
    Region = "Right Arm"
  ),
  CDQ009B = data.frame(
    Variable = c("PHQ-9", "PHQ-9:AD", "PHQ-9:CAD", "Marginal Effect"),
    Mean = c(1.0247276, 0.9991800, 0.9499295, 0.9849607),
    Lower = c(0.9711631, 0.9472787, 0.9016084, 0.9349283),
    Upper = c(1.081246, 1.053925, 1.000840, 1.037671),
    Region = "Right Chest"
  ),
  CDQ009C = data.frame(
    Variable = c("PHQ-9", "PHQ-9:AD", "PHQ-9:CAD", "Marginal Effect"),
    Mean = c(0.9819354, 0.9557197, 1.0852925, 1.0335917),
    Lower = c(0.9060085, 0.8731374, 0.9966001, 0.9546784),
    Upper = c(1.064225, 1.046113, 1.181878, 1.119028),
    Region = "Neck"
  ),
  CDQ009D = data.frame(
    Variable = c("PHQ-9", "PHQ-9:AD", "PHQ-9:CAD", "Marginal Effect"),
    Mean = c(1.1490472, 0.9822228, 0.8466896, 1.0068549),
    Lower = c(1.0635823, 0.9100638, 0.7828356, 0.9332860),
    Upper = c(1.2413798, 1.0601032, 0.9157519, 1.0862231),
    Region = "Upper Sternum"
  ),
  CDQ009E = data.frame(
    Variable = c("PHQ-9", "PHQ-9:AD", "PHQ-9:CAD", "Marginal Effect"),
    Mean = c(1.1297676, 1.0612211, 0.8784347, 1.0384291),
    Lower = c(1.0207360, 0.9584606, 0.7936171, 0.9384208),
    Upper = c(1.2504455, 1.1749990, 0.9723172, 1.1490955),
    Region = "Lower Sternum"
  ),
  CDQ009F = data.frame(
    Variable = c("PHQ-9", "PHQ-9:AD", "PHQ-9:CAD", "Marginal Effect"),
    Mean = c(0.9578713, 1.0295484, 1.0503839, 1.0020105),
    Lower = c(0.9191950, 0.9867414, 1.0052536, 0.9629802),
    Upper = c(0.998175, 1.074212, 1.097540, 1.042623),
    Region = "Left Chest"
  ),
  CDQ009G = data.frame(
    Variable = c("PHQ-9", "PHQ-9:AD", "PHQ-9:CAD", "Marginal Effect"),
    Mean = c(1.031356, 1.006010, 1.035412, 1.060844),
    Lower = c(0.8940819, 0.8703600, 0.8976527, 0.9205464),
    Upper = c(1.189706, 1.162802, 1.194312, 1.222523),
    Region = "Left Arm"
  ),
  CDQ009H = data.frame(
    Variable = c("PHQ-9", "PHQ-9:AD", "PHQ-9:CAD", "Marginal Effect"),
    Mean = c(1.0325968, 0.9883898, 1.0318884, 1.0546200),
    Lower = c(0.9234180, 0.8870933, 0.9179728, 0.9524612),
    Upper = c(1.154684, 1.101253, 1.159940, 1.167736),
    Region = "Epigastric"
  )
)

# Combine all models into a single data frame
all_data <- bind_rows(model_data)

# First, let's assume the data from both the original and new models is in two separate data frames
# Create a new data frame for the second set of results (new models)

new_data <- data.frame(
  Region = rep(c("Right Arm", "Right Chest", "Neck", "Upper Sternum", "Lower Sternum", "Left Chest", "Left Arm", "Epigastric"), each = 4),
  Variable = rep(c("PHQ-9", "PHQ-9:AD", "PHQ-9:CAD", "Marginal Effect"), times = 8),
  Mean = c(1.0257, 1.0177, 1.0166, 1.0048, 0.9750, 1.0039, 1.0181, 0.9570, 1.0105, 0.9611, 0.9647, 1.0804,
           1.0093, 1.1640, 0.9738, 0.8376, 1.0369, 1.1232, 1.0556, 0.8850, 1.0066, 0.9601, 1.0390, 1.0503,
           1.0444, 1.0167, 1.0211, 1.0285, 1.0377, 1.0099, 1.0073, 1.0370),
  Lower = c(0.9674, 0.9589, 0.9444, 0.9454, 0.9278, 0.9543, 0.9668, 0.9102, 0.9419, 0.8963, 0.8912, 1.0035,
               0.9423, 1.0845, 0.9087, 0.7795, 0.9435, 1.0212, 0.9599, 0.8046, 0.9679, 0.9201, 0.9955, 1.0022,
               0.9214, 0.8970, 0.8981, 0.9074, 0.9514, 0.9152, 0.9149, 0.9363),
  Upper = c(1.0875, 1.0801, 1.0942, 1.0679, 1.0246, 1.0560, 1.0723, 1.0061, 1.0842, 1.0306, 1.0443, 1.1631,
               1.0811, 1.2493, 1.0435, 0.9000, 1.1395, 1.2353, 1.1608, 0.9735, 1.0469, 1.0017, 1.0845, 1.1007,
               1.1838, 1.1525, 1.1611, 1.1658, 1.1318, 1.1143, 1.0872, 1.1485),
  Model = "Interactions-Only"
)

# Assuming all_data is the previous data frame that holds the original models
all_data$Model <- "Full Model"

# Combine both datasets
combined_data <- rbind(all_data, new_data)

ggplot(combined_data, aes(x = Mean, y = reorder(Variable, desc(Variable)), xmin = Lower, xmax = Upper, shape = Model, color = Model)) +
  geom_point(size = 3, position = position_dodge(width = 0.4)) +  # Adjust position for the points
  geom_errorbarh(height = 0.2, position = position_dodge(width = 0.4)) +  # Adjust position for the error bars
  facet_wrap(~Region, nrow = 2, ncol = 4, scales = "free_y") +
  theme_minimal() +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  geom_text(aes(label = round(Mean, 2)), vjust = -0.7, size = 3, color = "black", position = position_dodge(width = 0.4)) +  # Adjust text position
  scale_shape_manual(values = c(16, 15)) +  # Set point shapes: 16 for circle (original), 15 for square (new)
  scale_color_manual(values = c("#FF6666", "#9999FF")) +  # Set colors: blue for original, red for new
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 10, angle = 45, hjust = 1, vjust = 0, 
                               face = ifelse(combined_data$Variable == "Marginal Effect", "bold", "plain")),  # Bold Marginal Effect
    axis.text.x = element_text(size = 8),
    panel.spacing = unit(1.5, "lines"),
    strip.placement = "outside",
    panel.background = element_rect(fill = "white"),  # Background color for panels
    legend.title = element_blank(),  # Remove the legend title
    legend.text = element_text(size = 12),  # Increase the size of legend text
    legend.key.size = unit(1.5, "lines"),  # Increase the size of the legend key
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, lineheight = 0.9),  # Title formatting
    plot.margin = margin(t = 10, r = 10, b = 40, l = 10)  # Adjust the bottom margin (b = 40)
  ) +
  labs(
    x = "Odds-Ratio (95% CI)", 
    y = "Predictor",
    title = "Depression Score (PHQ-9) Odds-ratio Forest-Plot\nfor Logistic Regression on Chest Pain Location"
  )
