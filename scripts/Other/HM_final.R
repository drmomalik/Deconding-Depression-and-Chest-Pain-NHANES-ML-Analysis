
human_body_data <- hbd



## Plot for the full model 
# Ensure 'value' column is numeric (if not already)
human_body_data$value <- as.numeric(human_body_data$value)

human_body_data$value <- c(0, 1.15, 1.13, 0, 0.96, 0, 0, 0)


# Plot with a custom gradient
ggplot(human_body_data) +
  geom_sf(aes(fill = value), alpha = 0.9) +  # Add alpha for transparency
  scale_fill_gradient(
    low = "white",          # White for 0
    high = "red",           # Gradient for non-zero values
    na.value = "white"      # Ensure missing values are also white
  ) +
  geom_sf_text(
    aes(label = ifelse(value == 0, "", round(value, 2))),  # Show OR values for non-zero
    size = 7, color = "black"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),  # Remove gridlines for a cleaner map
    legend.title = element_text(size = 14),  # Adjust legend title size
    legend.text = element_text(size = 12),   # Adjust legend text size
    legend.key.size = unit(1.5, "cm")        # Adjust size of legend keys
  ) +
  labs(
    fill = "Odds Ratio"
  ) +
  theme(
    legend.position = "bottom"  # Move legend to bottom
  )


## Plot for the simple model 
# Ensure 'value' column is numeric (if not already)
human_body_data$value <- as.numeric(human_body_data$value)

human_body_data$value <- c(0, 1.16, 1.12, 0, 0, 0, 0, 0)


# Plot with a custom gradient
ggplot(human_body_data) +
  geom_sf(aes(fill = value), alpha = 0.9) +  # Add alpha for transparency
  scale_fill_gradient(
    low = "white",          # White for 0
    high = "red",           # Gradient for non-zero values
    na.value = "white"      # Ensure missing values are also white
  ) +
  geom_sf_text(
    aes(label = ifelse(value == 0, "", round(value, 2))),  # Show OR values for non-zero
    size = 7, color = "black"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),  # Remove gridlines for a cleaner map
    legend.title = element_text(size = 14),  # Adjust legend title size
    legend.text = element_text(size = 12),   # Adjust legend text size
    legend.key.size = unit(1.5, "cm")        # Adjust size of legend keys
  ) +
  labs(
    fill = "Odds Ratio"
  ) +
  theme(
    legend.position = "bottom"  # Move legend to bottom
  )

