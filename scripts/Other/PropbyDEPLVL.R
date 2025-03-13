library(ggplot2)
library(dplyr)
library(tidyr)

# Create a list of variables to plot
cdq_vars <- paste0("CDQ009", LETTERS[1:8])

# Function to calculate proportions for a given stratification variable
plot_proportions <- function(data, strat_var, title_prefix) {
  data %>%
    select(all_of(cdq_vars), all_of(strat_var)) %>%
    pivot_longer(cols = all_of(cdq_vars), 
                 names_to = "Variable", 
                 values_to = "Response") %>%
    filter(!is.na(Response)) %>%
    group_by(Variable, !!sym(strat_var), Response) %>%
    summarize(Count = n(), .groups = "drop") %>%
    group_by(Variable, !!sym(strat_var)) %>%
    mutate(Proportion = Count / sum(Count)) %>%
    ungroup() %>%
    ggplot(aes_string(x = "Response", y = "Proportion", fill = "as.factor(Response)")) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_grid(paste0("Variable ~ ", strat_var), scales = "free_y") +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(
      title = paste(title_prefix, "Proportions"),
      x = "Response",
      y = "Proportion",
      fill = "Response"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

# Generate plots for DEPR_LVL and DEPR_BIN
plot_depr_lvl <- plot_proportions(subsetd_2, "DEPR_LVL", "DEPR_LVL Stratified")
plot_depr_bin <- plot_proportions(subsetd_2, "DEPR_BIN", "DEPR_BIN Stratified")

# Print the plots
print(plot_depr_lvl)
print(plot_depr_bin)

