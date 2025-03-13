### Model of DEPR_TOT and TOT_REG 


imputation <- complete(imputations, 1)

# Define survey design
svy <- svydesign(
  id = ~SDMVPSU,
  weights = ~MEC15YR,
  strata = ~SDMVSTRA,
  data = imputation,
  nest = TRUE
)

# Subset to exclude rows where all response variables are 0 or NA
svy <- subset(svy, !apply(
  svy$variables[, responses], 
  1, 
  function(x) all(x == 0 | is.na(x))
))

# Fit logistic regression model for the current response
mod1 <- svyglm(TOT_REG ~ DEPR_TOT, design = svy, family = gaussian)
mod <- svyglm(TOT_REG ~ DEPR_TOT + AGE_BIN + RIAGENDR + RIDRETH1 + HIQ011 + I(BMXBMI^2) + 
                INC3 + DMDEDUC2 + PAQMV + I(ALQ130^2) + DIDTOT + DUQTOT + 
                SMQ040 * SMQ020 + DMDBORNT + BPQ020 + BPQ080 + SDDSRVYR + 
                DEPR_TOT * MEDDEP + DEPR_TOT * CADTOT + CDQ008 + CDQ010, design = svy, family = gaussian)
mod2 <- svyglm(TOT_REG ~ DEPR_TOT + AGE_BIN + RIAGENDR + RIDRETH1 + HIQ011 + I(BMXBMI^2) + 
                        INC3 + DMDEDUC2 + PAQMV + I(ALQ130^2) + DIDTOT + DUQTOT + 
                        SMQ040 * SMQ020 + DMDBORNT + BPQ020 + BPQ080 + SDDSRVYR + 
                       CDQ008 + CDQ010, design = svy, family = gaussian)

summary(mod)


x <- svy$variables$DEPR_TOT
y <- svy$variables$TOT_REG

# Scatter plot with explicit quadratic fit
ggplot(svy$variables, aes(x, y)) +
  geom_jitter(alpha = 0.6, color = "blue", size = 2, width = 0.3, height = 0.3) +  # Jittered points
  stat_smooth(
    method = "lm",
    formula = y ~ x,  
    se = TRUE, color = "red", size = 1
  ) +
  labs(
    title = "Depression Score vs Total Regions",
    x = "Depression Total Score",
    y = "Total Regions of Chest Pain"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  # Center the title
  )
