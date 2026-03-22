rm(list = ls())

plotTheme <- function() {
  theme(
    plot.title = element_text(size = 14, family = "sans", face = "plain", hjust = 0),
    plot.subtitle=element_text(size = 11, family = "sans", hjust = 0),
    plot.caption=element_text(size = 10, family = "sans", face = "italic", hjust = 0), 
    axis.title.x = element_text(size = 10, family = "sans", face = "plain", hjust = 1, vjust = -0.5),
    axis.title.y = element_text(size = 10, family = "sans", face = "plain", hjust = 1, vjust = 1),
    axis.text = element_text(size = 10, family = "sans", face = "plain"),
    panel.background = element_blank(),
    panel.grid.minor = element_line(colour = "gray"),
    panel.grid.major = element_line(colour = "gray"),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 10, family = "sans"),
    legend.text = element_text(size = 9, family = "sans"),
    axis.line = element_blank()
  )
}

library(nhanesA)
library(dplyr)
library(ggplot2)

## 1. Load NHANES data (2013–2014)
demo <- nhanes('DEMO_H')     # demographics
bmx  <- nhanes('BMX_H')      # body measurements

nhanes_data <- demo %>%
  select(SEQN, RIDAGEYR, RIAGENDR, RIDRETH1) %>%   # keep race/ethnicity
  inner_join(
    bmx %>% select(SEQN, BMXBMI, BMXWAIST, BMXARMC, BMXLEG), 
    by = "SEQN"
  ) %>%
  filter(!is.na(BMXBMI)) %>%                       # ensure BMI is observed
  mutate(
    gender = ifelse(RIAGENDR == 1, "Male", "Female"),
    race   = factor(RIDRETH1)                      # categorical
  )

nhanes_data <- nhanes_data[complete.cases(nhanes_data),]

Y <- nhanes_data$BMXBMI
N <- nrow(nhanes_data)

## 2. True population mean (NHANES sample mean)
true_mean <- mean(Y)

## 3. Make labeling informative: more likely for young participants
set.seed(123)
xi_true <- plogis(3 - 0.05 * nhanes_data$RIDAGEYR)  # in (0,1)
R <- rbinom(N, 1, xi_true)

n_lab <- sum(R)
cat("Number labeled:", n_lab, "\n")

## 4. Prediction model on labeled subset
fit <- lm(BMXBMI ~ RIDAGEYR + RIAGENDR + BMXWAIST + BMXARMC + BMXLEG + race, 
          data = nhanes_data[R == 1, ])
Yhat <- predict(fit, newdata = nhanes_data)

## 5. Estimates + CIs
z <- qnorm(0.975)

# Classic
classic <- mean(Y[R == 1])
se_classic <- sqrt(var(Y[R == 1]) / n_lab)

# HT (Equation 2.10)
ht <- mean(R * Y / xi_true)
var_ht <- (1/N^2) * sum(R * (1 - xi_true) / xi_true^2 * (Y - ht)^2)
se_ht <- sqrt(var_ht)

# Hájek (Equation 2.12)
w <- R / xi_true
w_sum <- sum(w)
hajek <- sum(w * Y) / w_sum
var_hajek <- sum((1 - xi_true) / xi_true^2 * R * (Y - hajek)^2) / (w_sum^2)
se_hajek <- sqrt(var_hajek)

# PPI unweighted
delta_unw <- mean(Yhat[R == 1] - Y[R == 1])
ppi_unw <- mean(Yhat) - delta_unw
s2_bias_unw <- mean(((Yhat[R == 1] - Y[R == 1]) - delta_unw)^2)
s2_yhat <- mean((Yhat - mean(Yhat))^2)
se_ppi_unw <- sqrt(s2_bias_unw / n_lab + s2_yhat / N)

# PPI weighted (Hájek-type variance for delta)
w_lab <- (R / xi_true)[R == 1]
xi_lab <- xi_true[R == 1]
e_lab <- (Yhat - Y)[R == 1]
delta_w <- sum(w_lab * e_lab) / sum(w_lab)
ppi_w <- mean(Yhat) - delta_w

# Variance of delta_w using linearization (Equation 2.12)
var_delta_w <- sum((1 - xi_lab) / xi_lab^2 * R[R==1] * (e_lab - delta_w)^2) / (sum(w_lab))^2
se_ppi_w <- sqrt(var_delta_w + s2_yhat / N)

## 6. Combine into results table
nhanes_results <- tibble(
  Estimator = c("Classic", "HT", "Hájek", "PPI_unweighted", "PPI_weighted"),
  Estimate = c(classic, ht, hajek, ppi_unw, ppi_w),
  SE = c(se_classic, se_ht, se_hajek, se_ppi_unw, se_ppi_w)
) %>%
  mutate(
    CI_lower = Estimate - z * SE,
    CI_upper = Estimate + z * SE,
    Bias = Estimate - true_mean,
    Width = CI_upper - CI_lower
  )

print(nhanes_results)

## 7. Plot stacked CIs
(ci_plot_nhanes <- ggplot(nhanes_results,
                          aes(y = Estimator,
                              xmin = CI_lower, xmax = CI_upper,
                              color = Estimator)) +
    geom_errorbarh(height = 0.2, size = 2, alpha = 0.5) +
    geom_point(aes(x = Estimate), size = 3) +
    geom_vline(xintercept = true_mean, linetype = "dashed", size = 1, color = "orange") +
    labs(x = "Estimated Mean BMI",
         y = NULL,
         caption = paste("Population mean BMI =", round(true_mean, 2))) +
    theme_bw() +
    theme(legend.position = "none") + plotTheme())

ggsave("art/nhanes_intervals_estimated.pdf", ci_plot_nhanes, width = 8, height = 6)

## 8. Summary statistics
summary_df <- nhanes_results %>%
  summarise(
    Estimator = Estimator,
    Bias = Bias,
    Width = Width
  )

print(summary_df)
