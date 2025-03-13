library(ggplot2)
library(splines)

data <- complete(imputations, 1)
complete_data <- complete(imputations, 1)
data <- subset(complete_data, !apply(
  complete_data[, responses], 1,
  function(x) all(x == 0 | is.na(x))
))
data <- data[,c(predictors, responses)]



# Plot splines with different `df`
ggplot(data, aes(DEPR_TOT, CDQ009A)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ ns(x, df = 3), color = "blue", se = FALSE) +
  geom_smooth(method = "lm", formula = y ~ ns(x, df = 5), color = "red", se = FALSE) +
  geom_smooth(method = "lm", formula = y ~ ns(x, df = 7), color = "green", se = FALSE) +
  labs(title = "Comparing Different Degrees of Freedom for Splines",
       subtitle = "Blue: df=3, Red: df=5, Green: df=7")
