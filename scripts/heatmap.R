library(here)
library(ggplot2)
library(tidyverse)
library(sf)
library(viridis)
library(lwgeom)



## Workflow
# Import the an SVG file into Inkscape and convert to dxf file
# Open QGIS, go to layer > add layer > add vector layer and select the dxf file
# create a new shapefile, use polygon to fill, add two fields: 
# one character field to name the region and one decimal or integer field for value/attributes
# Create all polygons then select the shapefile layer, right click, export
# safe as shapefile, then follow instruction below to import into R for plotting

## Will need to change the values in the dataframe to the actual observed values 
## being used for the heatmap (ex. odds ratio, proportion, proensity score)


# Load the shapefile (replace with your file path)
shapefile_path <- here("sf/body_part.shp")
hbd <- st_read(shapefile_path)

# > human_body_data
# body_parts value
# 1          neck     1
# 2 upper_sternum     1
# 3 lower_sternum     1
# 4   right_chest     1
# 5    left_chest     1
# 6     right_arm     1
# 7      left_arm     1
# 8    epigastric     1



human_body_data <- hbd

# No dep
human_body_data$value <- c(1.6, 5.6, 1.6, 5.6, 5.6, 1.6, 1.6, 1.6)
human_body_data$significant <- c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)

# Dep 
human_body_data <- hbd
human_body_data$value <- c(6.6, 2.2, 7, 1.6, 2, 4.6, 6.3, 1.9)
human_body_data$significant <- c(FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE)


# Ensure 'value' column is numeric (if not already)
human_body_data$value <- as.numeric(human_body_data$value)
set.seed(123)  # Optional: Set a seed for reproducibility
human_body_data$value <- runif(nrow(human_body_data), min = 0, max = 1)


# # Plot heatmap with ggplot2
# ggplot(data = human_body_data) +
#   geom_sf(aes(fill = value)) +   # Map the value to the fill color
#   scale_fill_viridis_c(option = "C") +  # Use a color scale (viridis is often a good choice)
#   theme_minimal() +                # Clean theme
#   labs(title = "Heatmap of Body Parts", 
#        fill = "Value") +           # Add legend title
#   theme(legend.position = "bottom")   # Adjust legend position




ggplot(human_body_data) +
  geom_sf(aes(fill = value), alpha = 0.9) +
  geom_sf_text(
    aes(label = sprintf("+%.1f%%", value)),  # Format value with 1 decimal place
    size = 5,                                # Adjust text size
    color = "black"                          # Improve readability against the heatmap
  ) +
  scale_fill_gradient(low = "beige", high = "red") +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.key.size = unit(1.5, "cm"),  # Increase legend key size
    legend.text = element_text(size = 14),  # Increase legend text size
    legend.title = element_text(size = 16)  # Increase legend title size
  ) +
  labs(fill = "Increase in\n probability\n of CAD (+%)") +
  geom_sf(data = human_body_data %>% filter(significant == TRUE), 
          color = "black", linewidth = 2,  alpha = 0.3, fill = NA)  # Thicker border with transparency


###  Text applied to only certain regions 
ggplot(human_body_data) +
  geom_sf(aes(fill = value), alpha = 0.9, color=NA) +
  geom_sf_text(
    aes(label = case_when(
      body_parts %in% c("epigastric", "upper_sternum") ~ sprintf("+%.1f%%", value),
      TRUE ~ ""  # Leave blank for other regions
    )), 
    size = 7, 
    color = "black"
  ) +
  scale_fill_gradient(low = "beige", high = "red") +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "none",  # Hide the legend
  ) +
  labs(fill = "Increase in\n probability\n of CAD (+%)")



### 


## plotting propotions in body heat map 
library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)

# Assuming `data` is your data frame containing the CDQ009* variables and DEPR_LVL
# `human_body_data` is your sf object with the geometry and body parts

# List of variables corresponding to body parts
body_parts <- c("right_arm", "right_chest", "neck", "upper_sternum", 
                "lower_sternum", "left_chest", "left_arm", "left_chest")
variables <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", 
               "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")

### Calculate proportions for each level of DEPR_LVL (example with just 1st level)
proportions <- df_lim %>%
  group_by(DEPR_LVL) %>%
  summarize(across(all_of(variables), ~ mean(. == 1, na.rm = TRUE), .names = "prop_{col}")) %>%
  pivot_longer(cols = starts_with("prop_"), 
               names_to = "variable", 
               values_to = "proportion") %>%
  mutate(variable = gsub("prop_", "", variable))

# Join proportions with `human_body_data`
human_body_data_1 <- human_body_data
human_body_data_1$value <- c(0.0907, 0.610, 0.180,
                             0.173, 0.363, 0.0659,
                             0.0563, 0.0124)


# Plot without the labels 
ggplot(human_body_data_1) +
geom_sf(aes(fill = value), alpha = 0.9) +  # Add alpha for transparency
scale_fill_gradient(low = "blue", high = "red") +
theme_minimal() +
theme(
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  axis.title = element_blank()) +
labs(fill = "Propotion")


## Basic log regression to build multiple heatmaps 
# List of variables corresponding to body parts
body_parts <- c("right_arm", "right_chest", "neck", "upper_sternum", 
                "lower_sternum", "left_chest", "left_arm", "epigastric")
variables <- c("CDQ009A", "CDQ009B", "CDQ009C", "CDQ009D", 
               "CDQ009E", "CDQ009F", "CDQ009G", "CDQ009H")

model <- glm(CDQ009A~-DEPR_LVL, data=wd, family=binomial, weights = MEC15YR)
summary(model)


