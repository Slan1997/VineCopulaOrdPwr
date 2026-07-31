############################################################
## Application 1:
## Simple Composition Example – Ranking by Clinical Severity
##
## Ordinal Composite Endpoint:
# defined as the most severe event occurring within 30 days of 
# emergency department admission and consists of six ordered categories, 
## W = 0,...,5
## 0: no event
## 1: ACS/PCI/CABG
## 2: emergent dialysis
## 3: intubation/mechanical ventilation
## 4: mechanical cardiac support
## 5: death
############################################################


source_local <- function(file) {
  if (file.exists(file)) {
    source(file)
  } else if (file.exists(file.path("R", file))) {
    source(file.path("R", file))
  } else {
    stop("Cannot find file: ", file)
  }
}

source_local("03_functions_vco_pmf_power.R")

library(pacman)
p_load(
  dplyr,
  tidyr,
  readr,
  stringr,
  magrittr,
  rvinecopulib,
  VineCopula # for BiCopTau2Par()
)

############################################################
# 1. Marginal component probabilities from STRATIFY
############################################################
p_component_control <- c(
  acs_pci_cabg = 0.059108527,
  emergent_dialysis = 0.020348837,
  intubation = 0.030038760,
  mechanical_cardiac_support = 0.005813953,
  death = 0.041666667
)

# Empirical ordinal composite distribution
p_ordinal_empirical <- c(
  0.875968992,
  0.046511628,
  0.013565891,
  0.016472868,
  0.005813953,
  0.041666667
)

############################################################
# 2. Sample size calculations
############################################################

get_cat_probs_from_control <- function(p_C, theta, A_Rate = 1, print = TRUE) {
  # p_C: category-specific probabilities in the control arm
  #      ordered from best to worst, or consistently with your PO model
  # theta: log odds ratio on the cumulative odds scale
  # A_Rate: allocation ratio A = n_C / n_T
  #         A_Rate = 1 for equal allocation
  
  # Basic checks
  if (any(p_C < 0)) {
    stop("All elements of p_C must be nonnegative.")
  }
  if (abs(sum(p_C) - 1) > 1e-6) {
    stop("p_C must sum to 1.")
  }
  if (length(p_C) < 2) {
    stop("p_C must contain at least two categories.")
  }
  
  L <- length(p_C)
  
  # Control-arm cumulative probabilities up to category L-1
  # gamma_C^ell = P(W <= ell | X = 0), ell = 1,...,L-1 in R indexing
  gamma_C <- cumsum(p_C)[1:(L - 1)]
  
  # Treatment-arm cumulative probabilities under proportional odds model
  gamma_T <- 1 / (1 + (1 / gamma_C - 1) * exp(-theta))
  
  # Convert treatment cumulative probabilities back to category probabilities
  p_T <- c(
    gamma_T[1],
    diff(gamma_T),
    1 - gamma_T[L - 1]
  )
  
  # Allocation-weighted average category-specific probabilities
  # A_Rate = n_C / n_T
  p_avg <- (p_T + A_Rate * p_C) / (A_Rate + 1)
  
  if (print) {
    cat(
      "p_C:   ", paste0(round(p_C, 4), collapse = ", "), "\n",
      "p_T:   ", paste0(round(p_T, 4), collapse = ", "), "\n",
      "p_avg: ", paste0(round(p_avg, 4), collapse = ", "), "\n"
    )
  }
  
  list(
    p_C = p_C,
    p_T = p_T,
    p_avg = p_avg,
    gamma_C = gamma_C,
    gamma_T = gamma_T
  )
}

calc_sample_size <- function(
    p_control = NULL,
    theta = log(1.5),
    power = .8,
    allocation = 1,
    p_avg = NULL
){
  
  if(is.null(p_avg)){
    p_avg <- get_cat_probs_from_control(
      p_C = p_control,
      theta = theta,
      A_Rate = allocation,
      print = FALSE
    )$p_avg
  }
  
  sum_cube <- sum(p_avg^3)
  
  n <- 
    3 * (allocation + 1)^2 *
    (qnorm(.975) + qnorm(power))^2 /
    (theta^2 * allocation * (1 - sum_cube))
  
  ceiling(n)
}

get_results <- function(
    p_control,
    theta = log(1.5),
    power = .8
){
  
  n <- calc_sample_size(
    p_control = p_control,
    theta = theta,
    power = power
  )
  
  list(
    sample_size = n
  )
}

############################################################
# 3. Baseline methods
############################################################

results <- tibble(
  Method = c(
    "WH-Empirical",
    "CD-0",
    "D1-Gaussian-LC",
    "D1-Gaussian-MC",
    "Naive"
  )
)

run_power_scenarios <- function(
    p_control = NULL,
    p_avg = NULL
){
  
  tibble(
    SS_power80_OR1.25 =
      calc_sample_size(
        p_control = p_control,
        p_avg = p_avg,
        theta = log(1.25),
        power = .8
      ),
    
    SS_power90_OR1.25 =
      calc_sample_size(
        p_control = p_control,
        p_avg = p_avg,
        theta = log(1.25),
        power = .9
      ),
    
    SS_power80_OR1.5 =
      calc_sample_size(
        p_control = p_control,
        p_avg = p_avg,
        theta = log(1.5),
        power = .8
      ),
    
    SS_power90_OR1.5 =
      calc_sample_size(
        p_control = p_control,
        p_avg = p_avg,
        theta = log(1.5),
        power = .9
      )
  )
}

# Gold-standard WH empirical
results[1,2:5] <-
  run_power_scenarios(p_ordinal_empirical)


# Naive uniform distribution
results[5,2:5] <-
  run_power_scenarios(
    p_avg = rep(1/6,6)
  )

############################################################
# 4. Vine copula ordinal distributions
############################################################

ordinal_from_vine <- function(
    corr_matrix,
    vine_type = "D",
    family = "Gaussian",
    components_order = 1:5
){
  
  idx <- components_order
  A <- corr_matrix
  
  
  family_id <- switch(
    family,
    Gaussian = 1,
    Clayton = 3,
    Gumbel = 4,
    Frank = 5,
    Joe = 6
  )
  
  
  if(vine_type == "D"){
    
    tau_list <- list(
      c(A[idx[1],idx[2]],
        A[idx[2],idx[3]],
        A[idx[3],idx[4]],
        A[idx[4],idx[5]]),
      
      c(A[idx[1],idx[3]],
        A[idx[2],idx[4]],
        A[idx[3],idx[5]]),
      
      c(A[idx[1],idx[4]],
        A[idx[2],idx[5]]),
      
      A[idx[1],idx[5]]
    )
    
  } else if(vine_type == "C"){
    
    tau_list <- list(
      A[idx[1],idx[2:5]],
      A[idx[2],idx[3:5]],
      A[idx[3],idx[4:5]],
      A[idx[4],idx[5]]
    )
    
  } else {
    
    stop("vine_type must be either 'D' or 'C'.")
    
  }
  
  
  
  family_list <- lapply(
    tau_list,
    function(x) rep(family_id,length(x))
  )
  
  
  # independence if tau = 0
  family_list <- Map(
    function(f,t){
      ifelse(t == 0, 0, f)
    },
    family_list,
    tau_list
  )
  
  
  param_list <- Map(
    function(f,t)
      BiCopTau2Par(
        family=f,
        tau=t
      ),
    family_list,
    tau_list
  )
  
  
  method_list <- lapply(
    family_list,
    function(x){
      ifelse(x==0,"Indep",family)
    }
  )
  
  
  setup <- VCO_Initialize(
    nComp = 5,
    prev = p_component_control[idx],
    vine_type = vine_type,
    method_list = method_list,
    rho_list = NULL,
    param_list = param_list
  )
  
  
  joint <- gen_jointpmf(setup)$jointpmf
  
  
  Ordinal_PMF(
    joint,
    method="order",
    severity_score_increasing=order(idx)
  )
}


############################################################
# 5. Correlation scenarios
############################################################


corr_level <- matrix(
  c(
    "1","Low","Very Low","Very Low","Very Low",
    "Low","1","Low","-Very Low","Low",
    "Very Low","Low","1","Low","Medium",
    "Very Low","-Very Low","Low","1","-Very Low",
    "Very Low","Low","Medium","-Very Low","1"
  ),
  5,
  byrow=TRUE,
  dimnames=list(
    names(p_component_control),
    names(p_component_control)
  )
)


convert_corr <- function(x){
  
  y <- matrix(NA,nrow(x),ncol(x))
  
  y[x=="Very Low"] <- .05
  y[x=="Low"] <- .2
  y[x=="Medium"] <- .4
  y[x=="-Very Low"] <- 0
  y[x=="1"] <- 1
  
  y
}


A0 <- convert_corr(corr_level)

A5 <- A0
A5[A5!=1] <- .4


############################################################
# 6. Vine results
############################################################


vine_results <- list(
  LC = ordinal_from_vine(A0),
  MC = ordinal_from_vine(A5)
)


results[2,2:5] <-
  run_power_scenarios(
    ordinal_from_vine(
      diag(5)
    )
  )


results[3,2:5] <-
  run_power_scenarios(vine_results$LC)


results[4,2:5] <-
  run_power_scenarios(vine_results$MC)

############################################################
# Final table
############################################################

results




############################################################
# PMFs derived from low and moderate correlation assumptions
# under D-vine Gaussian
############################################################

pmf_A0 <- ordinal_from_vine(
  corr_matrix = A0,
  vine_type = "D",
  family = "Gaussian",
  components_order = 1:5
)

pmf_A5 <- ordinal_from_vine(
  corr_matrix = A5,
  vine_type = "D",
  family = "Gaussian",
  components_order = 1:5
)

pmf_A0
pmf_A5
