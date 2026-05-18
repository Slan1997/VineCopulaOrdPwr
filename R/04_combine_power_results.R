# -----------------------------------------------------------------------------
# 04_combine_power_results.R
# -----------------------------------------------------------------------------
# Combine the large-n simulation outputs from 02_run_simulation_array.R and
# compute Whitehead power using:
#   1. the empirical ordinal composite PMF from the large simulation;
#   2. the naive uniform ordinal PMF;
#   3. nine VCO-derived ordinal composite PMFs.
#
# Main input:
#   - results0521_bign_1e5/sim_composite_mat_S2_<scenario>_0521.csv
#
# Main output:
#   - final_VCO_out_full0521_1e5.csv
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Source shared settings and functions
# -----------------------------------------------------------------------------

source_local <- function(file) {
  if (file.exists(file)) {
    source(file)
  } else if (file.exists(file.path("R", file))) {
    source(file.path("R", file))
  } else {
    stop("Cannot find file: ", file)
  }
}

source_local("00_settings_simulation.R")
source_local("01_functions_simulation.R")
source_local("03_functions_vco_pmf_power.R")


# -----------------------------------------------------------------------------
# Scenario-level power calculation
# -----------------------------------------------------------------------------

get_power_out_full <- function(scena,
                               power_n = c(500, 1000),
                               raw_dir = raw_output_dir) {
  if (is.na(scena) || scena < 1 || scena > nrow(scenarios)) {
    stop("scena must be an integer between 1 and ", nrow(scenarios), ".")
  }

  scenario_row <- scenarios[scena, ]

  beta_X <- beta_X_list[[scenario_row$beta_setup]]
  p_ct <- p_ct_list[[scenario_row$p_ct_setup]]
  p_trt <- get_p_trt(beta_X = beta_X, p_ct = p_ct)

  input_file <- file.path(
    raw_dir,
    paste0("sim_composite_mat_S2_", scena, "_0521.csv")
  )

  if (!file.exists(input_file)) {
    stop("Simulation output file not found: ", input_file)
  }

  sim_composite_mat <- readr::read_csv(input_file, show_col_types = FALSE)

  empirical_values_big_n <- sim_composite_mat %>%
    dplyr::filter(Comp_BinY != 1) %>%
    dplyr::group_by(Comp_BinY, n) %>%
    dplyr::summarise(
      dplyr::across(c(est_coef_Xs, dplyr::starts_with("p")), mean),
      .groups = "drop"
    )

  # The large simulated sample size is used only to estimate the underlying
  # treatment effect and ordinal composite PMF. Power is evaluated at power_n.
  empirical_values <- dplyr::bind_rows(
    lapply(
      power_n,
      function(n_eval) {
        empirical_values_big_n %>%
          dplyr::rename(bign = n) %>%
          dplyr::mutate(n = n_eval)
      }
    )
  )

  power_nonVCO <- empirical_values %>%
    dplyr::mutate(
      theta_R = abs(est_coef_Xs),
      sum_cube = apply(dplyr::select(., paste0("p", 0:n_bin_comp))^3, 1, sum),
      mu_beta1 = theta_R * sqrt(n^3 * (1 - sum_cube) / 12 / (n + 1)^2) - stats::qnorm(0.975),
      power_wh_empirical_dist = stats::pnorm(mu_beta1),
      sum_cube_naive = (1 / (n_bin_comp + 1))^3 * (n_bin_comp + 1),
      mu_beta1_naive = theta_R * sqrt(n^3 * (1 - sum_cube_naive) / 12 / (n + 1)^2) -
        stats::qnorm(0.975),
      power_wh_naive_dist = stats::pnorm(mu_beta1_naive)
    ) %>%
    dplyr::select(Comp_BinY, n, power_wh_empirical_dist, power_wh_naive_dist)

  VCO_out0 <- empirical_values %>%
    dplyr::transmute(
      Comp_BinY,
      n,
      mean_est_coef = est_coef_Xs
    )

  VCO_power <- matrix(
    NA_real_,
    nrow = nrow(VCO_out0),
    ncol = 9
  )
  colnames(VCO_power) <- paste0("power_wh_VCO", 1:9)

  VCO_out <- VCO_out0 %>%
    dplyr::bind_cols(as.data.frame(VCO_power)) %>%
    dplyr::arrange(n, Comp_BinY) %>%
    dplyr::mutate(chunk = rep(seq_len(nrow(.) / 3), each = 3)) %>%
    dplyr::relocate(chunk, .before = 1)

  for (VCsetting in 1:9) {
    SetupInfo_trt <- make_vco_setup(
      vco_setting = VCsetting,
      prev = p_trt,
      nComp = n_bin_comp
    )

    SetupInfo_ct <- make_vco_setup(
      vco_setting = VCsetting,
      prev = p_ct,
      nComp = n_bin_comp
    )

    jointpmf_trt <- gen_jointpmf(SetupInfo_trt)$jointpmf
    jointpmf_ct <- gen_jointpmf(SetupInfo_ct)$jointpmf
    jointpmf <- (jointpmf_trt + jointpmf_ct) / 2

    VC_ordinal_mat <- rbind(
      Ordinal_PMF(jointpmf, method = "sum"),
      Ordinal_PMF(jointpmf, method = "order", severity_score_increasing = n_bin_comp:1),
      Ordinal_PMF(jointpmf, method = "order", severity_score_increasing = 1:n_bin_comp)
    )

    for (ci in unique(VCO_out$chunk)) {
      n1 <- unique(VCO_out$n[VCO_out$chunk == ci])
      assume_effect_size <- VCO_out$mean_est_coef[VCO_out$chunk == ci]

      VCOpower <- whitehead_power(
        VC_ordinal_mat = VC_ordinal_mat,
        est_coef_Xs = assume_effect_size,
        n = n1
      )

      VCO_out[VCO_out$chunk == ci, paste0("power_wh_VCO", VCsetting)] <- VCOpower
    }
  }

  dplyr::left_join(
    power_nonVCO,
    VCO_out %>% dplyr::select(Comp_BinY, n, power_wh_VCO1:power_wh_VCO9),
    by = c("Comp_BinY", "n")
  ) %>%
    dplyr::mutate(scena = scena) %>%
    dplyr::relocate(scena, .before = 1)
}


# -----------------------------------------------------------------------------
# Run all scenarios and save combined output
# -----------------------------------------------------------------------------

final_VCO_out_full <- NULL

for (scena0 in seq_len(nrow(scenarios))) {
  message("Processing scenario ", scena0, " of ", nrow(scenarios), "...")

  final_VCO_out_full <- dplyr::bind_rows(
    final_VCO_out_full,
    get_power_out_full(scena = scena0 #, 
                       # raw_dir = "~/VineCopulaOrdPwr/simulated_data/"
                       )
  )
}

summary_output_file = "../simulated_data/final_VCO_out_full0521_1e5.csv"

readr::write_csv(final_VCO_out_full, summary_output_file)
message("Saved combined power results to: ", summary_output_file)
