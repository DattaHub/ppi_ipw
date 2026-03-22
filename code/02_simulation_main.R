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

library(ggplot2)
library(dplyr)

## Function to compute estimates and CIs (CORRECTED VARIANCES)
one_run <- function(Y_true, R, xi, Yhat, alpha = 0.05, bayes_alpha = 1) {
  N <- length(Y_true)
  n_labeled <- sum(R)
  z <- qnorm(1 - alpha/2)
  
  ## Classic (naïve)
  classic <- mean(Y_true[R == 1])
  se_classic <- sqrt(var(Y_true[R == 1]) / n_labeled)
  
  ## HT (CORRECTED - Equation 2.10)
  ht <- mean(R * Y_true / xi)
  var_ht <- (1/N^2) * sum(R * (1 - xi) / xi^2 * (Y_true - ht)^2)
  se_ht <- sqrt(var_ht)
  
  ## Hájek (CORRECTED - Equation 2.12)
  w <- R / xi
  w_sum <- sum(w)
  hajek <- sum(w * Y_true) / w_sum
  var_hajek <- sum((1 - xi) / xi^2 * R * (Y_true - hajek)^2) / (w_sum^2)
  se_hajek <- sqrt(var_hajek)
  
  ## PPI unweighted rectifier (Science paper simple formula)
  delta_unw <- mean(Yhat[R == 1] - Y_true[R == 1])
  ppi_unw <- mean(Yhat) - delta_unw
  s2_bias_unw <- mean(((Yhat[R == 1] - Y_true[R == 1]) - delta_unw)^2)
  s2_yhat <- mean((Yhat - mean(Yhat))^2)
  se_ppi_unw <- sqrt(s2_bias_unw / n_labeled + s2_yhat / N)
  
  ## PPI weighted rectifier (Hájek-style) (CORRECTED)
  w_lab <- (R / xi)[R == 1]
  xi_lab <- xi[R == 1]
  e_lab <- (Yhat - Y_true)[R == 1]
  delta_w <- sum(w_lab * e_lab) / sum(w_lab)
  ppi_w <- mean(Yhat) - delta_w
  
  # Variance using linearization (Equation 2.12)
  var_delta_w <- sum((1 - xi_lab) / xi_lab^2 * (e_lab - delta_w)^2) / (sum(w_lab))^2
  se_ppi_w <- sqrt(var_delta_w + s2_yhat / N)
  
  ## Bayes (Beta-Binomial)
  successes <- sum(Y_true[R == 1])
  failures <- n_labeled - successes
  bayes_est <- (successes + bayes_alpha) / (n_labeled + 2 * bayes_alpha)
  bayes_ci <- qbeta(c(alpha/2, 1 - alpha/2),
                    bayes_alpha + successes,
                    bayes_alpha + failures)
  
  tibble(
    Estimator = c("Classic", "HT", "Hájek", "PPI_unweighted", "PPI_weighted", "Bayes"),
    Estimate = c(classic, ht, hajek, ppi_unw, ppi_w, bayes_est),
    CI_lower = c(classic - z * se_classic,
                 ht - z * se_ht,
                 hajek - z * se_hajek,
                 ppi_unw - z * se_ppi_unw,
                 ppi_w - z * se_ppi_w,
                 bayes_ci[1]),
    CI_upper = c(classic + z * se_classic,
                 ht + z * se_ht,
                 hajek + z * se_hajek,
                 ppi_unw + z * se_ppi_unw,
                 ppi_w + z * se_ppi_w,
                 bayes_ci[2])
  )
}

## ---- Fixed super-population ----
set.seed(123)
N <- 500
X_pop <- rnorm(N)
Y_true_pop <- rbinom(N, 1, plogis(X_pop))
true_mean_global <- mean(Y_true_pop)

n_reps <- 200   # 200 replicates as mentioned in paper

results <- lapply(1:n_reps, function(rep) {
  xi_true <- plogis(0.5*X_pop)  # true inclusion probs for simulation
  R <- rbinom(N, 1, xi_true)  # sample using xi_true
  Yhat <- plogis(X_pop + rnorm(N, 0, 0.5))
  
  # --- Estimate xi empirically via logistic regression ---
  xi_hat <- fitted(glm(R ~ X_pop, family = binomial))
  
  one_run(Y_true_pop, R, xi_hat, Yhat) %>%
    mutate(Replicate = rep,
           TrueMean = true_mean_global)
}) %>% bind_rows()

## ---- Summary table ----
summary_stats <- results %>%
  group_by(Estimator) %>%
  summarise(
    Mean_Estimate = mean(Estimate),
    Bias = mean(Estimate) - true_mean_global,
    Mean_Width = mean(CI_upper - CI_lower),
    Coverage = mean(CI_lower <= true_mean_global & CI_upper >= true_mean_global),
    .groups = "drop"
  )

print(summary_stats)

## ---- Plot horizontal CIs for first 10 reps ----
## ---- Plot horizontal CIs for first 10 reps ----
plot_data <- results %>% filter(Replicate <= 10) %>% filter(Estimator != "Bayes")

(ci_plot <- ggplot(plot_data, aes(y = interaction(Replicate, Estimator),
                                  xmin = CI_lower, xmax = CI_upper,
                                  color = Estimator)) +
    geom_errorbarh(height = 0.3, size = 1.5, alpha = 0.5) +
    geom_point(aes(x = Estimate), size = 2) +
    geom_vline(xintercept = true_mean_global, linetype = "dashed", color = "orange") +
    scale_x_continuous(limits = c(0.40, 0.65)) +  # Match original range
    labs(x = "Estimated Mean",
         title = "Classic, HT, Hájek and PPI (unweighted & weighted)",  # Add title back
         caption = paste("True mean =", round(true_mean_global, 4))) +
    theme_minimal() +
    theme(legend.position = "top",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()))


ggsave("art/ci_plot_estimated_v2.pdf", ci_plot, width = 6, height = 6)

## ---- Width comparison plot ----
(width_plot <- ggplot(summary_stats %>% filter(Estimator != "Bayes"), 
                      aes(x = Estimator, y = Mean_Width, fill = Estimator)) +
   geom_col() +
   labs(y = "Mean CI Width", x = NULL,
        title = "Average 95% CI Width Comparison") +
   theme_minimal() +
   theme(legend.position = "none"))

ggsave("art/width_comparison_plot.pdf", width_plot, width = 10, height = 6)
