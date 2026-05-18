# -----------------------------------------------------------------------------
# 03_functions_vco_pmf_power.R
# -----------------------------------------------------------------------------
# Reusable functions for deriving ordinal composite endpoint PMFs from binary
# component prevalences and vine-copula dependence assumptions, and for computing
# Whitehead power using the derived ordinal PMFs.
#
# This file should contain reusable PMF/power functions only. Scenario settings,
# simulation constants, and copula parameter objects are defined in:
#   - 00_settings_simulation.R
#
# Main functions:
#   - VCO_Initialize(): initialize pair lists and PMF containers for C-/D-vines
#   - gen_jointpmf(): derive the joint PMF of binary components
#   - Ordinal_PMF(): map the binary joint PMF to ordinal composite PMFs
#   - whitehead_power(): compute Whitehead power for ordinal PMFs
#   - wh_power(): scalar Whitehead power calculation with allocation ratio
#   - make_vco_setup(): construct SetupInfo for one VCO assumption
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# VCO setup
# -----------------------------------------------------------------------------

VCO_Initialize <- function(nComp,
                           prev,
                           vine_type = "D",
                           method_list,
                           rho_list = NULL,
                           param_list = NULL) {
  if (!vine_type %in% c("C", "D")) {
    stop("vine_type must be either 'C' or 'D'.")
  }
  if (length(prev) != nComp) {
    stop("length(prev) must equal nComp.")
  }
  if (any(prev < 0 | prev > 1, na.rm = TRUE)) {
    stop("All prevalence values must be between 0 and 1.")
  }

  nTree <- nComp - 1
  qs <- 1 - prev

  if (!is.null(rho_list)) {
    eta_list <- vector("list", nTree)
  } else {
    eta_list <- NULL
  }

  if (vine_type == "D") {
    # Pair list for a D-vine.
    pair_list <- vector("list", nTree)
    names(pair_list) <- paste0("T", 2:nComp)

    for (pr in seq_len(nTree)) {
      if (pr == 1) {
        pair_list[[pr]] <- matrix(c(1:nTree, 2:nComp), byrow = TRUE, nrow = 2)
      } else {
        pair_list[[pr]] <- matrix(c(1:(2 * (nComp - pr))), nrow = 2)
      }
    }

    before_bar_joint_node <- vector("list", nTree)
    names(before_bar_joint_node) <- paste0("T", 2:nComp)

    for (l in seq_len(nTree)) {
      tree_idx <- l + 1
      temp <- matrix(NA, nrow = tree_idx, ncol = nComp - l)

      for (r in seq_len(tree_idx)) {
        temp[r, ] <- r:(nTree - l + r)
      }

      cont_num <- apply(temp, 2, function(x) paste0(x, collapse = ""))
      n_c <- length(cont_num)

      if (n_c > 2) {
        final_num <- c(cont_num[1], rep(cont_num[2:(n_c - 1)], each = 2), cont_num[n_c])
      } else {
        final_num <- cont_num
      }

      before_bar_joint_node[[l]] <- paste0("V", final_num)
    }

    after_bar_nodes <- vector("list", nTree - 1)
    names(after_bar_nodes) <- paste0("T", 2:(nComp - 1))

    joint_nodes_needed <- vector("list", nTree - 1)
    names(joint_nodes_needed) <- paste0("T", 3:nComp)

    for (l in seq_len(nTree - 1)) {
      tree_idx <- l + 1
      temp <- matrix(NA, nrow = l, ncol = nComp - tree_idx)

      for (r in seq_len(l)) {
        temp[r, ] <- (r + 1):(nTree - l + r)
      }

      cont_num <- apply(temp, 2, function(x) paste0(x, collapse = ""))
      joint_nodes_needed[[l]] <- paste0("V", cont_num)
      after_bar_nodes[[l]] <- paste0("V", rep(cont_num, each = 2))
    }
  } else {
    # Pair list for a C-vine.
    pair_list <- vector("list", nTree)
    names(pair_list) <- paste0("T", 2:nComp)

    for (pr in seq_len(nTree)) {
      pair_list[[pr]] <- matrix(
        c(rep(1, nComp - pr), 2:(nComp - pr + 1)),
        byrow = TRUE,
        nrow = 2
      )
    }

    before_bar_joint_node <- vector("list", nTree)
    names(before_bar_joint_node) <- paste0("T", 2:nComp)

    for (l in seq_len(nTree)) {
      before_bar_joint_node[[l]] <- paste0(
        "V",
        paste0(paste0(1:l, collapse = ""), (l + 1):nComp)
      )
    }

    after_bar_nodes <- vector("list", nTree - 1)
    names(after_bar_nodes) <- paste0("T", 2:(nComp - 1))

    joint_nodes_needed <- vector("list", nTree - 1)
    names(joint_nodes_needed) <- paste0("T", 3:nComp)

    for (l in seq_len(nTree - 1)) {
      after_bar_nodes[[l]] <- paste0("V", rep(paste0(1:l, collapse = ""), nComp - l))
      joint_nodes_needed[[l]] <- paste0("V", rep(paste0(1:l, collapse = ""), nComp - l - 1))
    }
  }

  # PMF container for single variables and all tree-level joint nodes.
  V_single <- matrix(c(qs, prev), nrow = 2, byrow = TRUE)
  colnames(V_single) <- paste0("V", 1:nComp)
  rownames(V_single) <- paste0("X", 0:1)

  list_joint_pmf <- list(V_single)

  for (ti in 2:(nTree + 1)) {
    Tidx <- paste0("T", ti)
    joints <- unique(before_bar_joint_node[[Tidx]])

    temp <- matrix(NA, nrow = 2^ti, ncol = length(joints))
    colnames(temp) <- joints
    rownames(temp) <- sort(
      paste0(
        "X",
        expand.grid(rep(list(0:1), ti)) |>
          apply(1, function(x) paste0(x, collapse = ""))
      )
    )

    list_joint_pmf <- c(list_joint_pmf, list(temp))
  }

  names(list_joint_pmf) <- paste0("T", 1:nComp)

  list(
    nComp = nComp,
    nTree = nTree,
    prev = prev,
    qs = qs,
    rho_list = rho_list,
    param_list = param_list,
    eta_list = eta_list,
    method_list = method_list,
    pair_list = pair_list,
    before_bar_joint_node = before_bar_joint_node,
    after_bar_nodes = after_bar_nodes,
    joint_nodes_needed = joint_nodes_needed,
    list_joint_pmf = list_joint_pmf
  )
}


# -----------------------------------------------------------------------------
# Bivariate copula and binary PMF helpers
# -----------------------------------------------------------------------------

## now have 5 different methods: Independent, Gaussian, Clayton, Frank, Gumbel.
## Note: Independent equals to Gaussian with 0 correlation.
## A good website to visualize copula contour: 
# https://copulatheque.shinyapps.io/copulas/ 

C_fun <- function(u1,
                  u2,
                  specific_method = NULL,
                  parC = NULL) {
  if (is.null(specific_method)) {
    specific_method <- "Indep"
  }

  if (specific_method == "Indep") {
    out <- u1 * u2
  } else {
    if (is.null(parC)) {
      stop("parC must be provided for non-independent copulas.")
    }

    if (specific_method == "Gaussian") {
      out <- VGAM::pbinorm(q1 = stats::qnorm(u1), q2 = stats::qnorm(u2), cov12 = parC)
    } else if (specific_method == "Clayton") {
      out <- max((u1^(-parC) + u2^(-parC) - 1)^(-1 / parC), 0)
    } else if (specific_method == "Frank") {
      out <- -1 / parC * log(
        1 + (exp(-parC * u1) - 1) * (exp(-parC * u2) - 1) / (exp(-parC) - 1)
      )
    } else if (specific_method == "Gumbel") {
      out <- exp(-(((-log(u1))^parC + (-log(u2))^parC)^(1 / parC)))
    } else {
      stop("Unknown copula method: ", specific_method)
    }
  }

  out
}

# C_fun(.8,.9, method='Indep')
# C_fun(.8,.9, method='Gaussian',parC=.1)
# C_fun(.8,.9, method='Clayton',parC=1)
# C_fun(.8,.9, method='Frank',parC=3)
# C_fun(.8,.9, method='Gumbel',parC=1.5)

## compare to VineCopula::BiCopCDF, works the same. Can use BiCopCDF when need more options.
# https://github.com/tnagler/VineCopula/blob/main/R/BiCopCDF.R
# library(VineCopula)
# BiCopCDF(u1=.8,u2=.9,family=0,par=NULL)
# BiCopCDF(u1=.8,u2=.9,family=1,par=.1)
# BiCopCDF(u1=.8,u2=.9,family=3,par=1)
# BiCopCDF(u1=.8,u2=.9,family=5,par=3)
# BiCopCDF(u1=.8,u2=.9,family=4,par=1.5)


# find_parC is now only valid for Gaussian copula
## if we know the correlation, we can find parameter for copula
find_parC <- function(p1,
                      p2,
                      rho,
                      specific_method) {
  q1 <- 1 - p1
  q2 <- 1 - p2
  target <- rho * sqrt(p1 * q1 * p2 * q2) + p1 * p2

  parC_list <- (1:1000) * 0.001
  try_get_target <- sapply(
    parC_list,
    function(x) C_fun(p1, p2, specific_method = specific_method, parC = x)
  )

  parC <- parC_list[which.min(abs(try_get_target - target))]
  if (parC == 1) parC <- 0.99
  parC
}

# pair copula for two binary components
BiCop_PMF <- function(p1,
                      p2,
                      parC,
                      specific_method,
                      print = FALSE) {
  parC0 <- parC

  if (print) message("Copula method: ", specific_method)
  if (print && specific_method != "Indep") message("Input parameter: ", parC0)

  q1 <- 1 - p1
  q2 <- 1 - p2

  biPMF <- -1
  ct <- 0

  while (!all(biPMF >= 0) && ct < 10) {
    ct <- ct + 1
    key <- C_fun(u1 = q1, u2 = q2, specific_method = specific_method, parC = parC)
    biPMF <- c(key, q1 - key, q2 - key, 1 - q1 - q2 + key)

    if (!all(biPMF >= 0)) {
      parC <- parC - 0.001
    }
  }
  ## sometimes the output pmf have some term less than 0
  if (!all(biPMF >= 0)) {
    warning("biPMF has at least one negative entry.")
  }
  if (print && !is.na(parC) && !identical(parC, parC0) && specific_method != "Indep") {
    message("Parameter used after adjustment: ", parC)
  }

  biPMF
}

# rho is only accepted for Gaussian copula
get_temp_bi_conditional_pmf <- function(prevalences,
                                        pair,
                                        rho = NULL,
                                        param1 = NULL, 
                                        method1  
                                        ) {
  param0 <- param1

  temp_bi_conditional_pmf <- matrix(NA, nrow = 4, ncol = ncol(pair))
  colnames(temp_bi_conditional_pmf) <- apply(
    pair,
    2,
    function(x) paste0("V", paste0(x, collapse = ""))
  )

  eta_vec <- rep(NA, ncol(pair))

  if (!is.null(rho) && length(rho) == 1) {
    rho <- rep(rho, ncol(pair))
  }
  if (is.null(param0)) {
    param1 <- rep(NA, ncol(pair))
  }
  if (!is.null(param0) && length(param0) == 1) {
    param1 <- rep(param1, ncol(pair))
  }
  if (length(method1) == 1) {
    method1 <- rep(method1, ncol(pair))
  }

  for (i in seq_len(ncol(pair))) {
    bi_idx <- paste0(c("V", pair[, i]), collapse = "")
    v1_idx <- pair[1, i]
    v2_idx <- pair[2, i]
    p1 <- as.numeric(prevalences[v1_idx])
    p2 <- as.numeric(prevalences[v2_idx])

    if (method1[i] == "Gaussian" && is.na(param1[i])) {
      eta_vec[i] <- find_parC(
        p1 = p1,
        p2 = p2,
        rho = rho[i],
        specific_method = method1[i]
      )
      param1[i] <- eta_vec[i]
    }

    temp_bi_conditional_pmf[, bi_idx] <- BiCop_PMF(
      p1 = p1,
      p2 = p2,
      parC = param1[i],
      specific_method = method1[i]
    )
  }

  list(
    temp_bi_conditional_pmf = temp_bi_conditional_pmf,
    eta_vec = eta_vec
  )
}


# -----------------------------------------------------------------------------
# Joint binary PMF and ordinal composite PMFs
# -----------------------------------------------------------------------------

gen_jointpmf <- function(SetupInfo, ...) {
  # Initialize
  nComp <- SetupInfo$nComp
  nTree <- SetupInfo$nTree
  prevalences <- SetupInfo$prev
  rho_list <- SetupInfo$rho_list
  param_list <- SetupInfo$param_list
  eta_list <- SetupInfo$eta_list
  method_list <- SetupInfo$method_list
  pair_list <- SetupInfo$pair_list
  before_bar_joint_node <- SetupInfo$before_bar_joint_node
  after_bar_nodes <- SetupInfo$after_bar_nodes
  joint_nodes_needed <- SetupInfo$joint_nodes_needed
  list_joint_pmf <- SetupInfo$list_joint_pmf

  for (idx in seq_len(nTree)) {
    if (!is.null(rho_list)) {
      if (is.list(rho_list)) {
        rho <- rho_list[[idx]]
      } else {
        ## rho_list can also be a vector if assume correlation are the same for each layer.
        rho <- rho_list[idx]
      }
    } else {
      rho <- NULL
    }

    if (is.list(method_list)) {
      method <- method_list[[idx]]
    } else {
      ### method_list can also be a vector if assume pair copula method are the same for each layer.
      method <- method_list[idx]
    }

    if (!is.null(param_list)) {
      if (is.list(param_list)) {
        param <- param_list[[idx]]
      } else {
        ## param_list can also be a vector if assume parameters are the same for each layer's bivariate copula models.
        param <- param_list[idx]
      }
    } else {
      param <- NULL
    }

    T_idx <- paste0("T", idx + 1) 
    T_idx_prev <- paste0("T", idx)
    pair <- pair_list[[T_idx]]

    if (idx == 1) {
      get_bi_result <- get_temp_bi_conditional_pmf(
        prevalences = prevalences,
        pair = pair,
        rho = rho,
        param1 = param,
        method1 = method
      )

      temp_bi_conditional_pmf <- get_bi_result[[1]]
      rownames(temp_bi_conditional_pmf) <- paste0(rep(c("X0", "X1"), each = 2), c("0", "1"))

      if (is.null(param) && all(method == "Gaussian")) {
        eta_list[[idx]] <- get_bi_result[[2]]
      }

      list_joint_pmf[[T_idx]] <- temp_bi_conditional_pmf
    } else {
      for (prev_row in seq_len(nrow(prevalences))) {
        get_bi_result <- get_temp_bi_conditional_pmf(
          prevalences = prevalences[prev_row, ],
          pair = pair,
          rho = rho,
          param1 = param,
          method1 = method
        )

        temp_bi_conditional_pmf <- get_bi_result[[1]]
        rownames(temp_bi_conditional_pmf) <- paste0(rep(c("X0", "X1"), each = 2), c("0", "1"))

        if (is.null(param) && all(method == "Gaussian")) {
          eta_list[[idx]] <- rbind(eta_list[[idx]], get_bi_result[[2]])
        }

        cond_Vs <- joint_nodes_needed[[T_idx]]
        cond_prob <- list_joint_pmf[[(idx - 1)]][prev_row, cond_Vs]
        joint_prob <- temp_bi_conditional_pmf
        joint_name <- colnames(list_joint_pmf[[T_idx]])
        cond_on <- gsub("X", "", rownames(prevalences)[prev_row])

        for (k in seq_along(joint_name)) {
          joint_var <- joint_name[k]
          cond_V <- gsub("V", "", cond_Vs[k])
          cond_V_digits <- strsplit(cond_V, split = "")[[1]]
          pos <- as.vector(
            sapply(
              cond_V_digits,
              function(x) grep(x, strsplit(joint_var, split = "")[[1]])
            )
          )

          target_rows <- sapply(
            rownames(list_joint_pmf[[T_idx]]),
            function(x) paste0(strsplit(x, split = "")[[1]][pos], collapse = "") == cond_on
          )

          list_joint_pmf[[T_idx]][target_rows, joint_var] <- joint_prob[, k] * cond_prob[k]
        }
      }
    }

    if (idx != nTree) {
      cjn <- before_bar_joint_node[[T_idx]]
      cn <- after_bar_nodes[[T_idx]]
      len <- length(cn)

      separate_cjn <- sapply(gsub("V", "", cjn), function(x) strsplit(x, split = "")[[1]])
      separate_cn <- sapply(gsub("V", "", cn), function(x) strsplit(x, split = "")[[1]])

      non_conded_digit <- rep(NA, len)
      non_conded_digit_pos <- rep(NA, len)

      for (s in seq_len(len)) {
        if (is.null(dim(separate_cn))) {
          sub <- separate_cn[s]
        } else {
          sub <- separate_cn[, s]
        }

        non_conded_digit[s] <- separate_cjn[!separate_cjn[, s] %in% sub, s]
        non_conded_digit_pos[s] <- which(separate_cjn[, s] == non_conded_digit[s])
      }

      prevalences <- matrix(NA, ncol = len, nrow = 2^idx)
      colnames(prevalences) <- paste(cjn, cn, sep = "|")
      rownames(prevalences) <- sort(
        paste0(
          "X",
          expand.grid(rep(list(0:1), idx)) |>
            apply(1, function(x) paste0(x, collapse = ""))
        )
      )

      for (r in rownames(prevalences)) {
        v_names <- matrix(NA, nrow = idx + 1, ncol = len)

        for (v in seq_len(len)) {
          v_names[non_conded_digit_pos[v], v] <- "1"
          v_names[-non_conded_digit_pos[v], v] <- strsplit(gsub("X", "", r), split = "")[[1]]
        }

        rjoint <- paste0("X", apply(v_names, 2, function(x) paste0(x, collapse = "")))
        conded_prob <- list_joint_pmf[[T_idx_prev]][r, cn] # after bar
        conding_prob <- apply(
          rbind(rjoint, cjn),
          2,
          function(x) list_joint_pmf[[T_idx]][x[1], x[2]]
        ) # before bar

        prevalences[r, ] <- conding_prob / conded_prob
      }
    } else {
      break
    }
  }

  list(
    list_joint_pmf = list_joint_pmf,
    eta_list = eta_list,
    jointpmf = list_joint_pmf[[length(list_joint_pmf)]]
  )
}

Ordinal_PMF <- function(jointpmf,
                        method,
                        severity_score_increasing = NULL) {
  nComp <- log2(nrow(jointpmf))
  entry <- gsub("X", "", rownames(jointpmf))

  if (method == "sum") {
    nc <- nchar(gsub("0", "", entry))
    Ordinal_sum_PMF <- matrix(NA, ncol = nComp + 1, nrow = 1)
    colnames(Ordinal_sum_PMF) <- paste0("X", 0:nComp)
    Ordinal_sum_PMF[1, ] <- sapply(0:nComp, function(x) sum(jointpmf[nc == x]))
    return(Ordinal_sum_PMF)
  }

  if (method != "order") {
    stop("method must be either 'sum' or 'order'.")
  }

  if (!is.vector(severity_score_increasing)) {
    stop("severity_score_increasing must be specified when method = 'order'.")
  }

  entry_new <- sapply(
    entry,
    function(x) paste0(strsplit(x, split = "")[[1]][severity_score_increasing], collapse = "")
  )

  jointpmf_cp <- jointpmf
  rownames(jointpmf_cp) <- entry_new

  Ordinal_ranking_PMF <- matrix(NA, ncol = nComp + 1, nrow = 1)
  colnames(Ordinal_ranking_PMF) <- paste0("X", 0:nComp)

  for (j in 0:nComp) {
    if (j == 0) {
      Ordinal_ranking_PMF[, "X0"] <- jointpmf_cp[paste0(rep(0, nComp), collapse = ""), ]
    } else {
      items <- entry_new[grep(paste0(1, paste0(rep(0, nComp - j), collapse = ""), "$"), entry_new)]
      Ordinal_ranking_PMF[, paste0("X", j)] <- sum(jointpmf_cp[items, ])
    }
  }

  Ordinal_ranking_PMF
}


# -----------------------------------------------------------------------------
# Whitehead power helpers
# -----------------------------------------------------------------------------

whitehead_power <- function(VC_ordinal_mat,
                            est_coef_Xs,
                            n) {
  theta_R <- sapply(est_coef_Xs, function(x) max(0.01, -x))
  sum_cube <- rowSums(VC_ordinal_mat^3)
  mu_beta1 <- theta_R * sqrt(n^3 * (1 - sum_cube) / 12 / (n + 1)^2) - stats::qnorm(0.975)
  stats::pnorm(mu_beta1)
}

wh_power <- function(theta_R,
                     n,
                     pr,
                     A) {
  sum_cube <- sum(pr^3)
  mu_beta1 <- theta_R * sqrt(A * n^3 * (1 - sum_cube) / 3 / (n + 1)^2 / (A + 1)^2) -
    stats::qnorm(0.975)
  stats::pnorm(mu_beta1)
}



make_vco_setup <- function(vco_setting,
                           prev,
                           nComp = 6) {
  if (vco_setting == 1) {
    return(VCO_Initialize(
      nComp = nComp,
      prev = prev,
      vine_type = "D",
      method_list = rep("Indep", nComp - 1),
      rho_list = NULL,
      param_list = NULL
    ))
  }

  if (vco_setting == 2) {
    return(VCO_Initialize(
      nComp = nComp,
      prev = prev,
      vine_type = "D",
      method_list = rep("Gaussian", nComp - 1),
      rho_list = NULL,
      param_list = Gaussian_param_HC_est
    ))
  }

  if (vco_setting == 3) {
    return(VCO_Initialize(
      nComp = nComp,
      prev = prev,
      vine_type = "C",
      method_list = rep("Gaussian", nComp - 1),
      rho_list = NULL,
      param_list = Gaussian_param_HC_est
    ))
  }

  if (vco_setting == 4) {
    return(VCO_Initialize(
      nComp = nComp,
      prev = prev,
      vine_type = "D",
      method_list = rep("Clayton", nComp - 1),
      rho_list = NULL,
      param_list = Clayton_param_HC_est
    ))
  }

  if (vco_setting == 5) {
    return(VCO_Initialize(
      nComp = nComp,
      prev = prev,
      vine_type = "D",
      method_list = rep("Frank", nComp - 1),
      rho_list = NULL,
      param_list = Frank_param_HC_est
    ))
  }

  if (vco_setting == 6) {
    return(VCO_Initialize(
      nComp = nComp,
      prev = prev,
      vine_type = "D",
      method_list = rep("Gaussian", nComp - 1),
      rho_list = NULL,
      param_list = Gaussian_param_LC_est
    ))
  }

  if (vco_setting == 7) {
    return(VCO_Initialize(
      nComp = nComp,
      prev = prev,
      vine_type = "C",
      method_list = rep("Gaussian", nComp - 1),
      rho_list = NULL,
      param_list = Gaussian_param_LC_est
    ))
  }

  if (vco_setting == 8) {
    return(VCO_Initialize(
      nComp = nComp,
      prev = prev,
      vine_type = "D",
      method_list = rep("Clayton", nComp - 1),
      rho_list = NULL,
      param_list = Clayton_param_LC_est
    ))
  }

  if (vco_setting == 9) {
    return(VCO_Initialize(
      nComp = nComp,
      prev = prev,
      vine_type = "D",
      method_list = rep("Frank", nComp - 1),
      rho_list = NULL,
      param_list = Frank_param_LC_est
    ))
  }

  stop("vco_setting must be an integer from 1 to 9.")
}
