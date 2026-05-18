# -----------------------------------------------------------------------------
# 00_settings_simulation.R
# -----------------------------------------------------------------------------
# Shared settings for the ordinal composite endpoint simulation study.
#
# This file defines the simulation constants, scenario grid, component
# prevalence settings, Kendall's tau settings, treatment effects, and copula
# parameter objects used across the simulation, post-processing, and plotting
# scripts.
#
# Intended use:
#   source("R/00_settings_simulation.R")
#
# Main scripts using this file:
#   - 02_run_simulation_array.R
#   - 04_combine_power_results.R
#   - 05_generate_tables_figures.R
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Packages
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(magrittr)
  library(ggplot2)
  library(rms)
  library(VineCopula)
  library(rvinecopulib)
  library(ggh4x)
})


# -----------------------------------------------------------------------------
# Global simulation constants
# -----------------------------------------------------------------------------

n_bin_comp <- 6
n_comp_outcomes <- 4

# Total sample size per simulated trial.
# In the main simulation study, each trial has n = 10,000 total subjects under
# equal allocation. This large n is used to estimate the true composite-outcome
# PMF and treatment effect, not directly as the trial sample size for power
# calculation.
ns <- 1e4

# Number of simulated trials per scenario.
nsim <- 1e5

# Output locations.
raw_output_dir <- file.path(getwd(), "results0521_bign_1e5")
summary_output_file <- file.path(getwd(), "final_VCO_out_full0521_1e5.csv")


# -----------------------------------------------------------------------------
# Component prevalence settings
# -----------------------------------------------------------------------------

# Control-arm component prevalence settings used in the simulation study.
# These correspond to the component-prevalence scenarios described in the paper:
#
#   high0    : Low component-prevalence setting
#              Components have low event probabilities, alternating 0.02 and 0.04.
#
#   even_out : Moderate component-prevalence setting
#              Components have moderate event probabilities, alternating 0.20 and 0.30.
#
#   high1    : High component-prevalence setting
#              Components have high event probabilities, alternating 0.70 and 0.80.
#
# Here, p_ct denotes the control-arm marginal prevalences for the six binary
# components used to construct the ordinal composite endpoint.

p_ct_list <- list(
  high0    = c(0.02, 0.04, 0.02, 0.04, 0.02, 0.04),
  even_out = c(0.20, 0.30, 0.20, 0.30, 0.20, 0.30),
  high1    = c(0.70, 0.80, 0.70, 0.80, 0.70, 0.80)
)

# # Paper-facing labels for the component-prevalence settings.
# p_ct_labels <- c(
#   high0 = "Low",
#   even_out = "Moderate",
#   high1 = "High"
# )


# -----------------------------------------------------------------------------
# Dependence settings
# -----------------------------------------------------------------------------

# True Kendall's tau values used to simulate data.
# High- and low-dependence settings are defined through tree-specific Kendall's
# tau values that decrease geometrically from Tree 1 to Tree 5 at a rate of 2/3.

tau_list <- list(
  high = round(0.6 * c(1, (2 / 3)^(1:4)), 2),
  low  = round(0.2 * c(1, (2 / 3)^(1:4)), 2)
)

# Assumed Kendall's tau values used for VCO-based power estimation.
# These are intentionally separated from tau_list because they represent the
# assumed dependence inputs used by the design-stage VCO methods, rather than
# the true dependence values used to generate the simulated data.

tau_list_est <- list(
  high = round(0.5 * c(1, (2 / 3)^(1:4)), 2),
  low  = round(0.3 * c(1, (2 / 3)^(1:4)), 2)
)

# # Paper-facing labels for dependence settings.
# corr_labels <- c(
#   high = "HC",
#   low = "LC"
# )


# -----------------------------------------------------------------------------
# Treatment-effect setting
# -----------------------------------------------------------------------------

# Component-specific treatment effects on the log-odds scale.
# The current simulation uses a common component-level effect, log(0.8), across
# all six binary components.

beta_X_list <- list(
  same = rep(log(0.8), n_bin_comp)
)


# -----------------------------------------------------------------------------
# Scenario grid
# -----------------------------------------------------------------------------

# True model labels:
#   C1-Gaussian : C-vine with Gaussian bicopulas
#   D1-Gaussian : D-vine with Gaussian bicopulas
#   D1-Clayton  : D-vine with Clayton bicopulas
#   D1-Frank    : D-vine with Frank bicopulas
#   D2          : D-vine with mixed bicopula families

scenarios <- expand.grid(
  p_ct_setup = c("high0", "even_out", "high1"),
  corr = c("high", "low"),
  beta_setup = "same",
  true_model = c(
    "C1-Gaussian",
    "D1-Gaussian",
    "D1-Clayton",
    "D1-Frank",
    "D2"
  ),
  stringsAsFactors = FALSE
) %>%
  mutate(scena_idx = row_number())


# -----------------------------------------------------------------------------
# Copula parameters for true simulation models
# -----------------------------------------------------------------------------

# Gaussian bicopulas: family = 1
Gaussian_param_HC_true <- VineCopula::BiCopTau2Par(
  family = rep(1, n_bin_comp - 1),
  tau = tau_list$high
)

Gaussian_param_LC_true <- VineCopula::BiCopTau2Par(
  family = rep(1, n_bin_comp - 1),
  tau = tau_list$low
)

# Clayton bicopulas: family = 3
Clayton_param_HC_true <- VineCopula::BiCopTau2Par(
  family = rep(3, n_bin_comp - 1),
  tau = tau_list$high
)

Clayton_param_LC_true <- VineCopula::BiCopTau2Par(
  family = rep(3, n_bin_comp - 1),
  tau = tau_list$low
)

# Frank bicopulas: family = 5
Frank_param_HC_true <- VineCopula::BiCopTau2Par(
  family = rep(5, n_bin_comp - 1),
  tau = tau_list$high
)

Frank_param_LC_true <- VineCopula::BiCopTau2Par(
  family = rep(5, n_bin_comp - 1),
  tau = tau_list$low
)

# Mixed D-vine bicopula setting, D2.
# Family coding follows VineCopula::BiCopTau2Par:
#   3 = Clayton, 4 = Gumbel, 5 = Frank.
Mix_param_HC_true <- list(
  VineCopula::BiCopTau2Par(family = c(3, 3, 5, 4, 4), tau = rep(tau_list$high[1], 5)),
  VineCopula::BiCopTau2Par(family = c(3, 5, 5, 4), tau = rep(tau_list$high[2], 4)),
  VineCopula::BiCopTau2Par(family = c(3, 5, 4), tau = rep(tau_list$high[3], 3)),
  VineCopula::BiCopTau2Par(family = c(3, 4), tau = rep(tau_list$high[4], 2)),
  VineCopula::BiCopTau2Par(family = c(5), tau = tau_list$high[5])
)

Mix_param_LC_true <- list(
  VineCopula::BiCopTau2Par(family = c(3, 3, 5, 4, 4), tau = rep(tau_list$low[1], 5)),
  VineCopula::BiCopTau2Par(family = c(3, 5, 5, 4), tau = rep(tau_list$low[2], 4)),
  VineCopula::BiCopTau2Par(family = c(3, 5, 4), tau = rep(tau_list$low[3], 3)),
  VineCopula::BiCopTau2Par(family = c(3, 4), tau = rep(tau_list$low[4], 2)),
  VineCopula::BiCopTau2Par(family = c(5), tau = tau_list$low[5])
)


# -----------------------------------------------------------------------------
# Copula parameters for assumed VCO power-estimation models
# -----------------------------------------------------------------------------

# These parameters are used for the nine VCO-based power-calculation settings
# in the post-simulation analysis.

Gaussian_param_HC_est <- VineCopula::BiCopTau2Par(
  family = rep(1, n_bin_comp - 1),
  tau = tau_list_est$high
)

Gaussian_param_LC_est <- VineCopula::BiCopTau2Par(
  family = rep(1, n_bin_comp - 1),
  tau = tau_list_est$low
)

Clayton_param_HC_est <- VineCopula::BiCopTau2Par(
  family = rep(3, n_bin_comp - 1),
  tau = tau_list_est$high
)

Clayton_param_LC_est <- VineCopula::BiCopTau2Par(
  family = rep(3, n_bin_comp - 1),
  tau = tau_list_est$low
)

Frank_param_HC_est <- VineCopula::BiCopTau2Par(
  family = rep(5, n_bin_comp - 1),
  tau = tau_list_est$high
)

Frank_param_LC_est <- VineCopula::BiCopTau2Par(
  family = rep(5, n_bin_comp - 1),
  tau = tau_list_est$low
)


# -----------------------------------------------------------------------------
# Method labels for VCO-based power estimation
# -----------------------------------------------------------------------------

vco_method_labels <- c(
  power_wh_empirical_dist = "WH-Empirical",
  power_wh_VCO1 = "WH-VCO-CD-0",
  power_wh_VCO2 = "WH-VCO-D1-Gaussian-HC",
  power_wh_VCO3 = "WH-VCO-C1-Gaussian-HC",
  power_wh_VCO4 = "WH-VCO-D1-Clayton-HC",
  power_wh_VCO5 = "WH-VCO-D1-Frank-HC",
  power_wh_VCO6 = "WH-VCO-D1-Gaussian-LC",
  power_wh_VCO7 = "WH-VCO-C1-Gaussian-LC",
  power_wh_VCO8 = "WH-VCO-D1-Clayton-LC",
  power_wh_VCO9 = "WH-VCO-D1-Frank-LC",
  power_wh_naive_dist = "WH-Naive"
)
