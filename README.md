# Prediction-Powered Inference with Inverse Probability Weighting

This repository contains R code to reproduce all simulations and data analyses in:

**Datta, J. and Polson, N. G. (2025).** *Prediction-Powered Inference with Inverse Probability Weighting.*  
Preprint available at [arXiv:2508.10149](https://arxiv.org/abs/2508.10149)

---

## Overview

This paper extends the Prediction-Powered Inference (PPI) framework to handle **informative labeling** (missing at random, or MAR) by incorporating inverse probability weighting (IPW). We show that PPI rectification admits a direct design-based interpretation as classical ratio estimators (Horvitz-Thompson and Hájek) applied to prediction residuals.

**Key Contributions:**
- Connects PPI to classical survey sampling estimators
- Provides variance estimation formulas using linearization (delta method)
- Demonstrates that weighted PPI achieves lower bias than unweighted PPI under informative labeling
- Real data analysis using NHANES 2013-2014 BMI data
- Extensive simulations comparing five estimation approaches

---

## Repository Structure
```
.
├── README.md
├── code/
│   ├── 01_nhanes_analysis.R              # Section 3.1: NHANES real data analysis
│   ├── 02_simulation_main.R              # Section 3.2: Main simulation (Table 1, Figure 1)
│   ├── 03_simulation_varying_plab.R      # Section 3.3: Simulation with varying p_lab
│   └── utils.R                           # Helper functions for estimators
└── figures/
    ├── nhanes_intervals_estimated.pdf    # Figure 3 (NHANES)
    ├── ci_plot_estimated.pdf             # Figure 1 (Simulation CIs)
    └── ci_plot_estimated_faceted.pdf     # Figure 2 (Varying p_lab)
```

---

## Files Description

### Main Analysis Scripts

1. **`01_nhanes_analysis.R`**
   - Real data application using NHANES 2013-2014
   - Estimates population mean BMI under informative labeling
   - Compares Classic, HT, Hájek, PPI (unweighted), and PPI (weighted) estimators
   - Generates Figure 3 in the paper
   - **Corrected variance formulas** for HT and Hájek estimators (Equations 2.10 and 2.12)

2. **`02_simulation_main.R`**
   - Main simulation study (N=500, 200 replicates)
   - Binary outcome with informative labeling
   - Generates Table 1 and Figure 1 in the paper
   - Demonstrates bias, CI width, and coverage properties

3. **`03_simulation_varying_plab.R`**
   - Large-scale simulation (N=10,000)
   - Varies labeled proportion p_lab ∈ {0.01, 0.02, 0.05}
   - Generates Table 2 and Figure 2 in the paper
   - Shows weighted PPI achieves consistently lower bias

---

## Key Methods Implemented

All scripts implement five estimators:

1. **Classic (naïve)**: Mean of labeled observations only
2. **Horvitz-Thompson (HT)**: IPW estimator with design-based variance
3. **Hájek**: Ratio estimator (normalized IPW)
4. **PPI (unweighted)**: Standard PPI rectification from Angelopoulos et al. (2023)
5. **PPI (weighted)**: Hájek-type weighted PPI rectification (our contribution)

**Variance Estimation:**
- Uses correct survey sampling variance formulas (Särndal et al., 1992)
- HT: Equation 2.10 with $(1-\hat{\xi}_i)/\hat{\xi}_i^2$ weighting
- Hájek: Equation 2.12 with linearization (delta method)
- Bootstrap procedure for estimated propensities

---

## Requirements
```r
# Required R packages
install.packages(c(
  "nhanesA",      # For NHANES data
  "dplyr",        # Data manipulation
  "ggplot2",      # Visualization
  "purrr"         # Functional programming
))
```

**R version:** >= 4.0.0

---

## Reproducing Results

### NHANES Analysis (Section 3.1)
```r
source("code/01_nhanes_analysis.R")
```
**Outputs:**
- Figure 3: `figures/nhanes_intervals_estimated.pdf`
- Console output with bias and CI widths

**Key Results:**
- Classic: bias −0.628, width 0.396
- HT: bias −0.111, width 0.229
- Hájek: bias −0.0277, width 0.229
- PPI (weighted): bias 0.0101, width 0.328

### Main Simulation (Section 3.2)
```r
source("code/02_simulation_main.R")
```
**Outputs:**
- Figure 1: `figures/ci_plot_estimated.pdf`
- Table 1 printed to console

**Key Results:**
- All IPW methods achieve 96-100% coverage
- HT and Hájek: width 0.092, coverage 96.5%
- Weighted PPI: width 0.092, coverage 97.5%

### Varying p_lab Simulation (Section 3.3)
```r
source("code/03_simulation_varying_plab.R")
```
**Outputs:**
- Figure 2: `figures/ci_plot_estimated_faceted.pdf`
- Table 2 printed to console

**Key Results:**
- Weighted PPI achieves lowest bias across all p_lab values
- At p_lab = 0.01: weighted PPI bias 0.003 vs unweighted 0.008

---

## Citation

If you use this code, please cite:
```bibtex
@article{datta2025ppi,
  title={Prediction-Powered Inference with Inverse Probability Weighting},
  author={Datta, Jyotishka and Polson, Nicholas G.},
  journal={arXiv preprint arXiv:2508.10149},
  year={2025}
}
```

---

## Related Work

This work builds on:
- **Angelopoulos et al. (2023)**: *Prediction-Powered Inference*, Science
  - Original PPI framework
  - Handles covariate shift in Section 4.2
- **Särndal, Swensson, and Wretman (1992)**: *Model Assisted Survey Sampling*
  - Classical reference for variance estimation
- **Royall and Pfeffermann (1982)**: *Balanced samples and robust Bayesian inference*
  - Bayesian interpretation of ratio estimators

---

## Authors

**Jyotishka Datta**  
Department of Statistics, Virginia Tech  
[jyotishka@vt.edu](mailto:jyotishka@vt.edu)

**Nicholas G. Polson**  
Booth School of Business, University of Chicago

---

## License

MIT License - see LICENSE file for details

---

## Notes

- **Variance Formulas**: This code uses the **corrected** variance estimators from Equations 2.10 and 2.12 in the paper, which properly account for the $(1-\xi_i)/\xi_i^2$ factor in the HT and Hájek variances.

- **NHANES Data**: The `nhanesA` package downloads data directly from CDC. No manual data download required.

- **Reproducibility**: All scripts use `set.seed()` for reproducible random number generation.

- **Computational Time**: 
  - NHANES analysis: ~30 seconds
  - Main simulation (200 reps): ~2 minutes
  - Varying p_lab (600 reps total): ~5 minutes

---

## Questions?

For questions about the code or methods, please open an issue on GitHub or contact the authors.
