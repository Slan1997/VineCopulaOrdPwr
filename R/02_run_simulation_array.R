# -----------------------------------------------------------------------------
# 02_run_simulation_array.R
# -----------------------------------------------------------------------------
# Run one array-task simulation scenario for the ordinal composite endpoint
# simulation study.
#
# Each SLURM array task corresponds to one row of the scenario grid defined in
# 00_settings_simulation.R. For that scenario, this script simulates nsim trials
# with n = 10,000 total subjects, estimates the empirical PMF of four composite
# outcomes, and fits the corresponding binary or proportional odds models.
#
# Required files:
#   - 00_settings_simulation.R
#   - 01_functions_simulation.R
#
# Main output:
#   - results0521_bign_1e5/sim_composite_mat_S2_<scenario>_0521.csv
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Shared settings and helper functions
# -----------------------------------------------------------------------------

# These checks make the script work when run either from the project root
# using Rscript R/02_run_simulation_array.R or from inside the R/ folder.
if (file.exists(file.path("R", "00_settings_simulation.R"))) {
  source(file.path("R", "00_settings_simulation.R"))
} else {
  source("00_settings_simulation.R")
}

if (file.exists(file.path("R", "01_functions_simulation.R"))) {
  source(file.path("R", "01_functions_simulation.R"))
} else {
  source("01_functions_simulation.R")
}


# -----------------------------------------------------------------------------
# Scenario index
# -----------------------------------------------------------------------------

# If running on SLURM, use SLURM_ARRAY_TASK_ID. If running interactively, a
# scenario index can be passed as the first command-line argument.
get_scenario_index <- function() {
  slurm_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
  args <- commandArgs(trailingOnly = TRUE)
  
  if (nzchar(slurm_id)) {
    return(as.integer(slurm_id))
  }
  
  if (length(args) >= 1) {
    return(as.integer(args[1]))
  }
  
  stop(
    "No scenario index provided. Use SLURM_ARRAY_TASK_ID or pass the ",
    "scenario index as the first command-line argument."
  )
}

scena <- get_scenario_index()

if (is.na(scena) || scena < 1 || scena > nrow(scenarios)) {
  stop("Scenario index must be an integer between 1 and ", nrow(scenarios), ".")
}

scenario_row <- scenarios[scena, ]
print(scenario_row)

beta_X <- beta_X_list[[scenario_row$beta_setup]]
p_ct <- p_ct_list[[scenario_row$p_ct_setup]]
tau <- tau_list[[scenario_row$corr]]
true_copula <- scenario_row$true_model
p_trt <- get_p_trt(beta_X = beta_X, p_ct = p_ct)


# -----------------------------------------------------------------------------
# Select true copula parameters for the current scenario
# -----------------------------------------------------------------------------

if (true_copula %in% c("C1-Gaussian", "D1-Gaussian")) {
  if (scenario_row$corr == "high") {
    parameters <- Gaussian_param_HC_true
  } else {
    parameters <- Gaussian_param_LC_true
  }
  
} else if (true_copula == "D1-Clayton") {
  if (scenario_row$corr == "high") {
    parameters <- Clayton_param_HC_true
  } else {
    parameters <- Clayton_param_LC_true
  }
  
} else if (true_copula == "D1-Frank") {
  if (scenario_row$corr == "high") {
    parameters <- Frank_param_HC_true
  } else {
    parameters <- Frank_param_LC_true
  }
  
} else if (true_copula == "D2") {
  if (scenario_row$corr == "high") {
    parameters <- Mix_param_HC_true
  } else {
    parameters <- Mix_param_LC_true
  }
  
} else {
  stop("Unknown true_copula value: ", true_copula)
}


# -----------------------------------------------------------------------------
# Simulate trials and estimate composite-outcome quantities
# -----------------------------------------------------------------------------

# For each scenario, simulate nsim trials. Each trial has n = ns total subjects
# under equal treatment allocation. For each simulated trial, construct four
# composite outcomes from the six binary components and estimate the empirical
# PMF and treatment-effect coefficient for each composite outcome.

nblock <- length(ns) * nsim

sim_composite_mat <- expand.grid(
  Comp_BinY = 1:n_comp_outcomes,
  sim = 1:nsim,
  n = ns
)

prob_mat <- matrix(
  NA_real_,
  nrow = nrow(sim_composite_mat),
  ncol = n_bin_comp + 1
)
colnames(prob_mat) <- paste0("p", 0:n_bin_comp)

sim_composite_mat <- sim_composite_mat %>%
  mutate(
    block_idx = rep(1:nblock, each = n_comp_outcomes),
    est_coef_Xs = NA_real_,
    p_val_Xs = NA_real_
  ) %>%
  bind_cols(as.data.frame(prob_mat))

head(sim_composite_mat)

for (bi in 1:nblock) {
  if (bi %% 500 == 0) print(bi)
  
  row_idx <- which(sim_composite_mat$block_idx == bi)
  
  n0 <- unique(sim_composite_mat$n[row_idx])
  n1 <- n0 / 2  # treatment:control = 1:1
  
  set.seed(bi)
  
  if (true_copula == "D1-Gaussian") {
    dt_trt <- VC_sim_bin(
      n_1arm = n1,
      parameters = parameters,
      prevalence = p_trt,
      x_input = 1,
      vine = "D-vine",
      bicop = "gaussian",
      num_bin_comp = n_bin_comp
    )
    
    dt_ct <- VC_sim_bin(
      n_1arm = n1,
      parameters = parameters,
      prevalence = p_ct,
      x_input = 0,
      vine = "D-vine",
      bicop = "gaussian",
      num_bin_comp = n_bin_comp
    )
    
  } else if (true_copula == "C1-Gaussian") {
    dt_trt <- VC_sim_bin(
      n_1arm = n1,
      parameters = parameters,
      prevalence = p_trt,
      x_input = 1,
      vine = "C-vine",
      bicop = "gaussian",
      num_bin_comp = n_bin_comp
    )
    
    dt_ct <- VC_sim_bin(
      n_1arm = n1,
      parameters = parameters,
      prevalence = p_ct,
      x_input = 0,
      vine = "C-vine",
      bicop = "gaussian",
      num_bin_comp = n_bin_comp
    )
    
  } else if (true_copula == "D1-Clayton") {
    dt_trt <- VC_sim_bin(
      n_1arm = n1,
      parameters = parameters,
      prevalence = p_trt,
      x_input = 1,
      vine = "D-vine",
      bicop = "clayton",
      num_bin_comp = n_bin_comp
    )
    
    dt_ct <- VC_sim_bin(
      n_1arm = n1,
      parameters = parameters,
      prevalence = p_ct,
      x_input = 0,
      vine = "D-vine",
      bicop = "clayton",
      num_bin_comp = n_bin_comp
    )
    
  } else if (true_copula == "D1-Frank") {
    dt_trt <- VC_sim_bin(
      n_1arm = n1,
      parameters = parameters,
      prevalence = p_trt,
      x_input = 1,
      vine = "D-vine",
      bicop = "frank",
      num_bin_comp = n_bin_comp
    )
    
    dt_ct <- VC_sim_bin(
      n_1arm = n1,
      parameters = parameters,
      prevalence = p_ct,
      x_input = 0,
      vine = "D-vine",
      bicop = "frank",
      num_bin_comp = n_bin_comp
    )
    
  } else if (true_copula == "D2") {
    dt_trt <- VC_sim_bin(
      n_1arm = n1,
      parameters = parameters,
      prevalence = p_trt,
      x_input = 1,
      vine = "D-vine",
      bicop = "mixed",
      num_bin_comp = n_bin_comp
    )
    
    dt_ct <- VC_sim_bin(
      n_1arm = n1,
      parameters = parameters,
      prevalence = p_ct,
      x_input = 0,
      vine = "D-vine",
      bicop = "mixed",
      num_bin_comp = n_bin_comp
    )
  }
  
  full_dt_yx <- rbind(dt_trt, dt_ct)
  full_dt <- full_dt_yx[, 1:n_bin_comp]
  
  # Four composite outcomes:
  #   comp_Y1: any component event vs none
  #   comp_Y2: total number of component events
  #   comp_Y3: most severe event, where component 1 is most severe
  #   comp_Y4: most severe event, where component 6 is most severe
  comp_Y1 <- as.numeric(rowSums(full_dt) > 0)
  comp_Y2 <- as.factor(rowSums(full_dt))
  comp_Y3 <- as.factor(apply(full_dt %*% diag(n_bin_comp:1), 1, max))
  comp_Y4 <- as.factor(apply(full_dt %*% diag(1:n_bin_comp), 1, max))
  
  # Empirical PMF of each composite outcome.
  pr_Y1 <- c(
    sapply(0:1, function(x) mean(comp_Y1 == x)),
    rep(NA_real_, n_bin_comp - 1)
  )
  pr_Y2 <- sapply(0:n_bin_comp, function(x) mean(comp_Y2 == x))
  pr_Y3 <- sapply(0:n_bin_comp, function(x) mean(comp_Y3 == x))
  pr_Y4 <- sapply(0:n_bin_comp, function(x) mean(comp_Y4 == x))
  
  sim_composite_mat[row_idx, colnames(prob_mat)] <- rbind(
    pr_Y1,
    pr_Y2,
    pr_Y3,
    pr_Y4
  )
  
  dt_X <- data.frame(X = full_dt_yx[, "X"])
  
  # Fit models for the four composite outcomes.
  # comp_Y1 is binary. comp_Y2 may be binary in rare-event settings, so use a
  # binary model if only two categories are observed; otherwise use ordinal lrm.
  simY <- comp_Y1
  out1 <- fit_model(
    simY = simY,
    dt_X = dt_X,
    outcome_type = "binary",
    Print = FALSE
  )
  
  simY <- comp_Y2
  if (length(levels(comp_Y2)) == 2) {
    out2 <- fit_model(
      simY = simY,
      dt_X = dt_X,
      outcome_type = "binary",
      Print = FALSE
    )
  } else {
    dd <- rms::datadist(data.frame(dt_X, simY))
    options(datadist = "dd")
    
    out2 <- fit_model(
      simY = simY,
      dt_X = dt_X,
      outcome_type = "ordinal",
      Print = FALSE
    )
  }
  
  simY <- comp_Y3
  dd <- rms::datadist(data.frame(dt_X, simY))
  options(datadist = "dd")
  
  out3 <- fit_model(
    simY = simY,
    dt_X = dt_X,
    outcome_type = "ordinal",
    Print = FALSE
  )
  
  simY <- comp_Y4
  dd <- rms::datadist(data.frame(dt_X, simY))
  options(datadist = "dd")
  
  out4 <- fit_model(
    simY = simY,
    dt_X = dt_X,
    outcome_type = "ordinal",
    Print = FALSE
  )
  
  sim_composite_mat[row_idx, "est_coef_Xs"] <- c(
    out1$est_beta_X,
    out2$est_beta_X,
    out3$est_beta_X,
    out4$est_beta_X
  )
  
  sim_composite_mat[row_idx, "p_val_Xs"] <- c(
    out1$p_value_X,
    out2$p_value_X,
    out3$p_value_X,
    out4$p_value_X
  )
}

head(sim_composite_mat)

dir.create(raw_output_dir, recursive = TRUE, showWarnings = FALSE)

output_file <- file.path(
  raw_output_dir,
  paste0("sim_composite_mat_S2_", scena, "_0521.csv")
)

readr::write_csv(sim_composite_mat, output_file)