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
library(purrr)

## ----- helper: calibrate intercept so mean(plogis(a + b*X)) = target -----
calibrate_intercept <- function(X, slope = 0.5, target = 0.02) {
  f <- function(a) mean(plogis(a + slope * X)) - target
  uniroot(f, interval = c(-15, 15))$root
}

## ----- estimators (CORRECTED VARIANCES) -----
one_run <- function(Y_true, R, xi, Yhat, alpha = 0.05, bayes_alpha = 1) {
  N <- length(Y_true)
  n_labeled <- sum(R)
  z <- qnorm(1 - alpha/2)
  
  ## Classic
  classic <- mean(Y_true[R == 1])
  se_classic <- sqrt(stats::var(Y_true[R == 1]) / n_labeled)
  
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
  
  ## PPI (unweighted rectifier)
  delta_unw <- mean(Yhat[R == 1] - Y_true[R == 1])
  ppi_unw <- mean(Yhat) - delta_unw
  s2_bias_unw <- mean(((Yhat[R == 1] - Y_true[R == 1]) - delta_unw)^2)
  s2_yhat <- mean((Yhat - mean(Yhat))^2)
  se_ppi_unw <- sqrt(s2_bias_unw / n_labeled + s2_yhat / N)
  
  ## PPI (Hájek-style weighted rectifier) (CORRECTED)
  w_lab <- (R / xi)[R == 1]
  xi_lab <- xi[R == 1]
  e_lab <- (Yhat - Y_true)[R == 1]
  delta_w <- sum(w_lab * e_lab) / sum(w_lab)
  ppi_w <- mean(Yhat) - delta_w
  
  # Variance using linearization (Equation 2.12)
  var_delta_w <- sum((1 - xi_lab) / xi_lab^2 * (e_lab - delta_w)^2) / (sum(w_lab))^2
  se_ppi_w <- sqrt(var_delta_w + s2_yhat / N)
  
  ## Bayes (Beta-Binomial) — unchanged
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

## ======= High N / tiny n_labeled simulation =======
set.seed(123)

N_pop   <- 10000                 # big population
X_pop   <- rnorm(N_pop)
Y_true_pop <- rbinom(N_pop, 1, plogis(X_pop))
true_mean_global <- mean(Y_true_pop)

## predictive model quality
Yhat <- plogis(X_pop + rnorm(N_pop, 0, 0.5))

## choose labeled fractions to explore: 1%, 2%, 5%
p_lab_grid <- c(0.01, 0.02, 0.05)
n_reps <- 200

results <- map_dfr(p_lab_grid, function(p_lab) {
  a_int <- calibrate_intercept(X_pop, slope = 0.5, target = p_lab)
  xi_true <- plogis(a_int + 0.5 * X_pop)          # true inclusion probs with small mean
  n_lab_expected <- round(mean(xi_true) * N_pop)
  
  cat(sprintf("\nRunning p_lab = %.2f (expected n ≈ %d)\n", p_lab, n_lab_expected))
  
  map_dfr(1:n_reps, function(rep) {
    R <- rbinom(N_pop, 1, xi_true)
    
    ## safety: if zero labels occur (rare at very small p_lab), resample
    tries <- 0
    while (sum(R) < 10 && tries < 10) {  # ensure at least some labels
      R <- rbinom(N_pop, 1, xi_true)
      tries <- tries + 1
    }
    
    ## estimate xi via logistic regression on full features (R observed for all units)
    xi_hat <- fitted(glm(R ~ X_pop, family = binomial))
    
    out <- one_run(Y_true_pop, R, xi_hat, Yhat)
    out %>%
      mutate(Replicate = rep,
             TrueMean = true_mean_global,
             p_lab = p_lab,
             n_lab = sum(R),
             ratio_unlab_to_lab = (N_pop - sum(R)) / sum(R))
  })
})

## ---- Summary table by labeled fraction ----
summary_stats <- results %>%
  group_by(p_lab, Estimator) %>%
  summarise(
    Mean_Estimate = mean(Estimate),
    Bias          = mean(Estimate) - first(TrueMean),
    Mean_Width    = mean(CI_upper - CI_lower),
    Coverage      = mean(CI_lower <= first(TrueMean) & CI_upper >= first(TrueMean)),
    Avg_n_lab     = mean(n_lab),
    Avg_ratio     = mean(ratio_unlab_to_lab),
    .groups = "drop"
  ) %>%
  arrange(p_lab, Estimator)

print(summary_stats)

## ---- CI plots for a single small-labeled regime (say p_lab = 0.02) ----
plot_data_single <- results %>% filter(p_lab == 0.02, Replicate <= 10, Estimator != "Bayes")

(ci_plot <- ggplot(plot_data_single,
                   aes(y = interaction(Replicate, Estimator),
                       xmin = CI_lower, xmax = CI_upper,
                       color = Estimator)) +
    geom_errorbarh(height = 0.3, size = 1.3, alpha = 0.6) +
    geom_point(aes(x = Estimate), size = 1.8) +
    geom_vline(xintercept = true_mean_global, linetype = "dashed", color = "orange") +
    labs(x = "Estimated Mean",
         caption = paste("True mean =", round(true_mean_global, 4))) +
    theme_minimal() +
    theme(legend.position = "top",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()))

ggsave("art/ci_plot_estimated_high_ratio_v2.pdf", ci_plot, width = 8, height = 5)

## ---- Faceted plot across p_lab values ----
(ci_plot_faceted <- ggplot(results %>% 
                             filter(Replicate <= 10, Estimator != "Bayes"),
                           aes(y = interaction(Replicate, Estimator),
                               xmin = CI_lower, xmax = CI_upper,
                               color = Estimator)) +
   geom_errorbarh(height = 0.3, size = 1.3, alpha = 0.6) +
   geom_point(aes(x = Estimate), size = 1.8) +
   geom_vline(xintercept = true_mean_global, 
              linetype = "dashed", color = "orange") +
   labs(x = "Estimated Mean",
        caption = paste("True mean =", round(true_mean_global, 4))) +
   theme_minimal() +
   theme(legend.position = "top",
         axis.text.y = element_blank(),
         axis.ticks.y = element_blank()) +
   facet_grid(. ~ p_lab, labeller = label_both)
)

ggsave("art/ci_plot_estimated_faceted_v2.pdf", ci_plot_faceted, width = 14, height = 5)

## ---- Width comparison plot ----
(width_plot <- ggplot(plot_data_single,
                      aes(x = Estimator, y = CI_upper - CI_lower, fill = Estimator)) +
   geom_col(width = 0.7) +
   labs(y = "CI Width", x = "Estimator",
        title = "Confidence Interval Widths (p_lab = 0.02)") +
   theme_minimal() +
   theme(legend.position = "none"))

ggsave("art/ci_widths_high_ratio.pdf", width_plot, width = 7, height = 4)
