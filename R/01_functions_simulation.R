# -----------------------------------------------------------------------------
# 01_functions_simulation.R
# -----------------------------------------------------------------------------
# Helper functions for simulating correlated binary components using vine copulas
# and fitting regression models to the resulting composite outcomes.
#
# Main functions used by the simulation script:
#   - get_p_trt(): compute treatment-arm component prevalences under log-ORs
#   - VC_sim_bin(): simulate correlated binary components for one trial arm
#   - fit_model(): estimate treatment effect and p-value for a composite outcome
#
# Internal helper functions:
#   - odds(): convert probabilities to odds
#   - setup_VC(): construct the vine copula object used by VC_sim_bin()
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Probability helpers
# -----------------------------------------------------------------------------

odds <- function(p) {
  if (any(p <= 0 | p >= 1, na.rm = TRUE)) {
    stop("All probabilities must be strictly between 0 and 1.")
  }
  p / (1 - p)
}

get_p_trt <- function(beta_X, p_ct) {
  if (length(beta_X) != length(p_ct)) {
    stop("beta_X and p_ct must have the same length.")
  }

  exp_beta <- exp(beta_X)
  exp_beta * odds(p_ct) / (exp_beta * odds(p_ct) + 1)
}


# -----------------------------------------------------------------------------
# Vine copula setup for simulating true data
# -----------------------------------------------------------------------------

# vine_type:
#   "C-vine" or "D-vine"
#
# bicop_type:
#   family names used by rvinecopulib::bicop_dist(), such as
#   "gaussian", "clayton", "frank", or "indep";
#   or "mixed" for the six-component D2 simulation setting.

setup_VC <- function(vine_type,
                     bicop_type,
                     parameters,
                     num_bin_comp = 6) {
  if (!vine_type %in% c("C-vine", "D-vine")) {
    stop("vine_type must be either 'C-vine' or 'D-vine'.")
  }
  if (num_bin_comp < 2) {
    stop("num_bin_comp must be at least 2.")
  }

  n_tree <- num_bin_comp - 1

  if (bicop_type != "mixed") {
    if (length(parameters) != n_tree) {
      stop("For a non-mixed bicopula model, length(parameters) must equal num_bin_comp - 1.")
    }

    pcs <- vector("list", n_tree)
    for (tree in seq_len(n_tree)) {
      bicop_temp <- rvinecopulib::bicop_dist(
        family = bicop_type,
        rotation = 0,
        parameters = parameters[tree],
        var_types = c("d", "d")
      )
      pcs[[tree]] <- rep(list(bicop_temp), num_bin_comp - tree)
    }

    if (vine_type == "D-vine") {
      mat <- rvinecopulib::dvine_structure(num_bin_comp:1)
    } else {
      mat <- rvinecopulib::cvine_structure(num_bin_comp:1) # the order does not matter since using all same bicops
    }
  } else {
    # Mixed bicopula model used for the D2 simulation setting.
    # This implementation is currently specific to the six-component D-vine.
    if (num_bin_comp != 6 || vine_type != "D-vine") {
      stop("The mixed bicopula setup is currently implemented only for a six-component D-vine.")
    }
    if (!is.list(parameters) || length(parameters) != 5) {
      stop("For bicop_type = 'mixed', parameters must be a list of length 5.")
    }

    mixed_families <- list(
      c("clayton", "clayton", "frank", "gumbel", "gumbel"),
      c("clayton", "frank", "frank", "gumbel"),
      c("clayton", "frank", "gumbel"),
      c("clayton", "gumbel"),
      c("frank")
    )

    pcs <- vector("list", length(mixed_families))
    for (tree in seq_along(mixed_families)) {
      if (length(parameters[[tree]]) != length(mixed_families[[tree]])) {
        stop("Each element of parameters must match the number of bicopulas in the corresponding tree.")
      }

      pcs[[tree]] <- Map(
        function(fam, par) {
          rvinecopulib::bicop_dist(
            family = fam,
            rotation = 0,
            parameters = par,
            var_types = c("d", "d")
          )
        },
        fam = mixed_families[[tree]],
        par = parameters[[tree]]
      )
    }

    mat <- rvinecopulib::dvine_structure(6:1)
  }

  rvinecopulib::vinecop_dist(
    pair_copulas = pcs,
    structure = mat,
    var_types = rep("d", num_bin_comp)
  )
}


# -----------------------------------------------------------------------------
# Simulation of correlated binary components
# -----------------------------------------------------------------------------

VC_sim_bin <- function(n_1arm, # sample size for one arm
                       parameters,
                       prevalence,
                       x_input, # 1: treatment; 0: control
                       vine,
                       bicop,
                       num_bin_comp = 6) {
  if (length(prevalence) != num_bin_comp) {
    stop("length(prevalence) must equal num_bin_comp.")
  }
  if (any(prevalence < 0 | prevalence > 1, na.rm = TRUE)) {
    stop("All prevalence values must be between 0 and 1.")
  }

  vc_model <- setup_VC(
    vine_type = vine,
    bicop_type = bicop,
    parameters = parameters,
    num_bin_comp = num_bin_comp
  )

  cdata <- rvinecopulib::rvinecop(n = n_1arm, vinecop = vc_model)
  
  bin_dt <- cdata
  for (j in seq_len(ncol(cdata))) {
    cutoff <- stats::quantile(cdata[, j], probs = prevalence[j], names = FALSE)
    bin_dt[cdata[, j] < cutoff, j] <- 1
    bin_dt[cdata[, j] >= cutoff, j] <- 0
  }

  cbind(bin_dt, X = rep(x_input, n_1arm))
}


# -----------------------------------------------------------------------------
# Model fitting
# -----------------------------------------------------------------------------

fit_model <- function(simY,
                      dt_X,
                      outcome_type,
                      Print = FALSE) {
  if (!"X" %in% names(dt_X)) {
    stop("dt_X must contain a column named 'X'.")
  }

  analysis_dt <- data.frame(dt_X, simY = simY)

  if (outcome_type == "binary") {
    if (Print) message("Fitting logistic regression.")

    fit <- stats::glm(
      simY ~ X,
      family = stats::binomial(),
      data = analysis_dt
    )

    coef_tab <- summary(fit)$coefficients
    est_beta_X <- coef_tab[2,1]
    p_X <- coef_tab[2,4]

  } else if (outcome_type == "ordinal") {
    if (Print) message("Fitting proportional odds model.")

    fit <- rms::lrm(
      simY ~ X,
      data = analysis_dt
    )

    est_beta_X <- as.numeric(fit$coefficients["X"])
    p_X <- anova(fit)["X", "P"]
  } else {
    stop("outcome_type must be one of: 'binary' or 'ordinal'.")
  }

  list(
    est_beta_X = est_beta_X,
    p_value_X = p_X
  )
}
