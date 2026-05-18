# -----------------------------------------------------------------------------
# 05_generate_tables_figures.R
# -----------------------------------------------------------------------------
# Generate manuscript-ready tables and figures from the combined power results.
#
# Input:
#   - final_VCO_out_full0521_1e5.csv
#     produced by 04_combine_power_results.R
#
# Output:
#   - tables/
#   - figures/
#
# Notes:
#   - True model names are kept as full family names, e.g., D1-Gaussian.
#   - Estimator/method labels are shortened, e.g., D1-Gaus-HC.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Source shared settings
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


# -----------------------------------------------------------------------------
# User-controlled output settings
# -----------------------------------------------------------------------------

# The manuscript tables/figures below use one evaluation sample size at a time.
# Change this to 1000 if needed.
target_n <- 500

table_dir <- file.path("../tables")
figure_dir <- file.path("../figures")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

summary_output_file = "../simulated_data/final_VCO_out_full0521_1e5.csv"

if (!file.exists(summary_output_file)) {
  stop(
    "Combined power summary file not found: ", summary_output_file,
    "\nRun 04_combine_power_results.R first."
  )
}

# -----------------------------------------------------------------------------
# Display labels
# -----------------------------------------------------------------------------

# Correlation labels
corr_labels <- c(
  high = "HC",
  low  = "LC"
)

# Component prevalence labels
p_ct_labels <- c(
  high0    = "Low",
  even_out = "Moderate",
  high1    = "High"
)

# Ordinal composite type labels
comp_labels <- c(
  `2` = "Summation",
  `3` = "Rank 1",
  `4` = "Rank 2"
)

# True model labels: keep full family names
true_model_labels <- c(
  "C1-Gaussian" = "C1-Gaussian",
  "D1-Gaussian" = "D1-Gaussian",
  "D1-Clayton"  = "D1-Clayton",
  "D1-Frank"    = "D1-Frank",
  "D2"          = "D2"
)

# Method labels: keep WH-Empirical, shorten copula family names for others
method_labels <- c(
  power_wh_empirical_dist = "WH-Empirical",
  power_wh_VCO2           = "D1-Gaus-HC",
  power_wh_VCO3           = "C1-Gaus-HC",
  power_wh_VCO4           = "D1-Clay-HC",
  power_wh_VCO5           = "D1-Fran-HC",
  power_wh_VCO6           = "D1-Gaus-LC",
  power_wh_VCO7           = "C1-Gaus-LC",
  power_wh_VCO8           = "D1-Clay-LC",
  power_wh_VCO9           = "D1-Fran-LC",
  power_wh_VCO1           = "CD-0",
  power_wh_naive_dist     = "Naive"
)

method_order_all <- unname(method_labels)


# -----------------------------------------------------------------------------
# Read and format combined power results
# -----------------------------------------------------------------------------

final_VCO_out_full <- readr::read_csv(
  summary_output_file,
  show_col_types = FALSE
) %>%
  dplyr::filter(n == target_n)

if (nrow(final_VCO_out_full) == 0) {
  stop("No rows found for target_n = ", target_n, ".")
}


# -----------------------------------------------------------------------------
# Join scenario labels and create manuscript-facing table
# -----------------------------------------------------------------------------

scenario_labels <- scenarios %>%
  dplyr::rename(scena = scena_idx) %>%
  dplyr::mutate(
    `True Correlation` = unname(corr_labels[corr]),
    `True Model` = unname(true_model_labels[true_model]),
    `True Settings` = paste0(`True Model`, "-", `True Correlation`),
    `Component Prev Type` = unname(p_ct_labels[p_ct_setup])
  ) %>%
  dplyr::select(
    scena,
    `True Settings`,
    `True Correlation`,
    `True Model`,
    `Component Prev Type`
  )

out_table <- final_VCO_out_full %>%
  dplyr::filter(Comp_BinY %in% 2:4) %>%
  dplyr::left_join(scenario_labels, by = "scena") %>%
  dplyr::mutate(
    `Ordinal Composite Type` = unname(comp_labels[as.character(Comp_BinY)])
  ) %>%
  dplyr::select(
    `Ordinal Composite Type`,
    `True Settings`,
    `True Correlation`,
    `True Model`,
    `Component Prev Type`,
    dplyr::all_of(names(method_labels))
  ) %>%
  dplyr::rename(!!!setNames(names(method_labels), method_labels)) %>%
  dplyr::arrange(
    factor(`Ordinal Composite Type`, levels = unname(comp_labels)),
    `Component Prev Type`,
    `True Settings`
  )

# Percent version for manuscript tables
out_table_percent <- out_table %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(method_order_all),
      ~ round(.x * 100, 1)
    )
  )

readr::write_csv(
  out_table,
  file.path(table_dir, paste0("power_table_all_methods_n", target_n, ".csv"))
)

readr::write_csv(
  out_table_percent,
  file.path(table_dir, paste0("power_table_all_methods_n", target_n, "_percent.csv"))
)


# -----------------------------------------------------------------------------
# Long-format table for plotting
# -----------------------------------------------------------------------------

out_long <- out_table %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(method_order_all),
    names_to = "Method",
    values_to = "Power"
  ) %>%
  dplyr::mutate(
    Method = factor(Method, levels = method_order_all),
    `Ordinal Composite Type` = factor(
      `Ordinal Composite Type`,
      levels = unname(comp_labels)
    ),
    `Component Prev Type` = factor(
      `Component Prev Type`,
      levels = c("Low", "Moderate", "High")
    ),
    `True Correlation` = factor(
      `True Correlation`,
      levels = c("HC", "LC")
    )
  )

# -----------------------------------------------------------------------------
# Helper functions for tables and plots
# -----------------------------------------------------------------------------

save_comparison_table <- function(data, comparison_name, comp_type) {
  clean_name <- comparison_name %>%
    stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
    stringr::str_to_lower()
  
  comp_name <- comp_type %>%
    stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
    stringr::str_to_lower()
  
  readr::write_csv(
    data,
    file.path(table_dir, paste0("out_table_", clean_name, "_", comp_name, ".csv"))
  )
}

make_power_plot <- function(data,
                            comp_type,
                            true_model,
                            methods_to_keep,
                            file_name,
                            y_limits = c(7, 57)) {
  plot_data <- data %>%
    dplyr::filter(
      `Ordinal Composite Type` == comp_type,
      `True Model` == true_model,
      Method %in% methods_to_keep
    ) %>%
    dplyr::mutate(
      Method = factor(Method, levels = methods_to_keep),
      `Component Prevalence Type` = "Component Prevalence Type",
      `Component Prev Type` = factor(
        `Component Prev Type`,
        levels = c("Low", "Moderate", "High")
      ),
      `True Correlation Label` = "True Correlation",
      `True Correlation` = factor(
        `True Correlation`,
        levels = c("HC", "LC"),
        labels = c("high", "low")
      )
    )
  
  if (nrow(plot_data) == 0) {
    warning("No data found for plot: ", file_name)
    return(invisible(NULL))
  }
  
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = Method, y = Power * 100)
  ) +
    ggplot2::geom_point(
      size = 1.2,
      stroke = 0.3
    ) +
    ggh4x::facet_nested(
      rows = ggplot2::vars(`True Correlation Label`, `True Correlation`),
      cols = ggplot2::vars(`Component Prevalence Type`, `Component Prev Type`)
    ) +
    ggplot2::labs(
      title = paste0("True Model: ", true_model),
      x = "Methods",
      y = "Power (%)"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1,
        size = 8
      ),
      axis.text.y = ggplot2::element_text(size = 9),
      axis.title.x = ggplot2::element_text(size = 10),
      axis.title.y = ggplot2::element_text(size = 10),
      
      # No borders around each panel
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(
        fill = "white",
        colour = NA
      ),
      
      # Keep only left and bottom axis lines
      axis.line.x.bottom = ggplot2::element_line(
        linewidth = 0.3,
        colour = "black"
      ),
      axis.line.y.left = ggplot2::element_line(
        linewidth = 0.3,
        colour = "black"
      ),
      axis.line.x.top = ggplot2::element_blank(),
      axis.line.y.right = ggplot2::element_blank(),
      
      panel.grid.major = ggplot2::element_line(
        linewidth = 0.25,
        colour = "grey85"
      ),
      panel.grid.minor = ggplot2::element_line(
        linewidth = 0.15,
        colour = "grey92"
      ),
      
      panel.spacing.x = grid::unit(0.3, "lines"),
      panel.spacing.y = grid::unit(0.3, "lines"),
      
      strip.background = ggplot2::element_rect(
        fill = "grey85",
        colour = NA
      ),
      strip.text = ggplot2::element_text(
        size = 9,
        colour = "grey20"
      ),
      
      strip.switch.pad.grid = grid::unit(0, "pt"),
      strip.switch.pad.wrap = grid::unit(0, "pt"),
      strip.clip = "off",
      ggh4x.facet.nestline = ggplot2::element_blank(),
      
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = 11,
        face = "bold"
      ),
      
      legend.position = "none"
    )
  
  if (!is.null(y_limits)) {
    p <- p + ggplot2::coord_cartesian(ylim = y_limits)
  }
  
  ggplot2::ggsave(
    filename = file.path(figure_dir, file_name),
    plot = p,
    width = 5.2,
    height = 5.0,
    dpi = 300
  )
  
  invisible(p)
}

# -----------------------------------------------------------------------------
# Round 1: Compare C-vine vs D-vine Gaussian specifications
# -----------------------------------------------------------------------------

methods_round1 <- c(
  "WH-Empirical",
  "C1-Gaus-HC",
  "D1-Gaus-HC",
  "C1-Gaus-LC",
  "D1-Gaus-LC",
  "CD-0",
  "Naive"
)

table_round1 <- out_table %>%
  dplyr::filter(`True Model` %in% c("C1-Gaussian", "D1-Gaussian")) %>%
  dplyr::select(
    `Ordinal Composite Type`,
    `True Settings`,
    `True Correlation`,
    `True Model`,
    `Component Prev Type`,
    dplyr::all_of(methods_round1)
  )

for (comp_type in unname(comp_labels)) {
  for (true_model_i in c("C1-Gaussian", "D1-Gaussian")) {
    
    file_stub <- paste0(
      "round1_",
      stringr::str_to_lower(stringr::str_replace_all(true_model_i, "-", "_")),
      "_",
      stringr::str_to_lower(stringr::str_replace_all(comp_type, " ", "_"))
    )
    
    readr::write_csv(
      table_round1 %>%
        dplyr::filter(
          `Ordinal Composite Type` == comp_type,
          `True Model` == true_model_i
        ),
      file.path(table_dir, paste0(file_stub, ".csv"))
    )
    
    make_power_plot(
      data = out_long,
      comp_type = comp_type,
      true_model = true_model_i,
      methods_to_keep = methods_round1,
      file_name = paste0(file_stub, ".png")
    )
  }
}


######## Round 2
# -----------------------------------------------------------------------------
# Comparison 2: D1-Clayton true model
# -----------------------------------------------------------------------------

methods_round2_clayton <- c(
  "WH-Empirical",
  "D1-Clay-HC",
  "D1-Gaus-HC",
  "D1-Clay-LC",
  "D1-Gaus-LC",
  "CD-0",
  "Naive"
)

table_round2_clayton <- out_table %>%
  dplyr::filter(`True Model` == "D1-Clayton") %>%
  dplyr::select(
    `Ordinal Composite Type`,
    `True Settings`,
    `True Correlation`,
    `True Model`,
    `Component Prev Type`,
    dplyr::all_of(methods_round2_clayton)
  )

for (comp_type in unname(comp_labels)) {
  
  file_stub <- paste0(
    "round2_d1_clayton_",
    stringr::str_to_lower(stringr::str_replace_all(comp_type, " ", "_"))
  )
  
  readr::write_csv(
    table_round2_clayton %>%
      dplyr::filter(`Ordinal Composite Type` == comp_type),
    file.path(table_dir, paste0(file_stub, ".csv"))
  )
  
  make_power_plot(
    data = out_long,
    comp_type = comp_type,
    true_model = "D1-Clayton",
    methods_to_keep = methods_round2_clayton,
    file_name = paste0(file_stub, ".png")
  )
}

# -----------------------------------------------------------------------------
# Comparison 3: D1-Frank true model
# -----------------------------------------------------------------------------

methods_round2_frank <- c(
  "WH-Empirical",
  "D1-Fran-HC",
  "D1-Gaus-HC",
  "D1-Fran-LC",
  "D1-Gaus-LC",
  "CD-0",
  "Naive"
)

table_round2_frank <- out_table %>%
  dplyr::filter(`True Model` == "D1-Frank") %>%
  dplyr::select(
    `Ordinal Composite Type`,
    `True Settings`,
    `True Correlation`,
    `True Model`,
    `Component Prev Type`,
    dplyr::all_of(methods_round2_frank)
  )

for (comp_type in unname(comp_labels)) {
  
  file_stub <- paste0(
    "round2_d1_frank_",
    stringr::str_to_lower(stringr::str_replace_all(comp_type, " ", "_"))
  )
  
  readr::write_csv(
    table_round2_frank %>%
      dplyr::filter(`Ordinal Composite Type` == comp_type),
    file.path(table_dir, paste0(file_stub, ".csv"))
  )
  
  make_power_plot(
    data = out_long,
    comp_type = comp_type,
    true_model = "D1-Frank",
    methods_to_keep = methods_round2_frank,
    file_name = paste0(file_stub, ".png")
  )
}


# -----------------------------------------------------------------------------
# Comparison 4: D2 mixed-family true model
# -----------------------------------------------------------------------------
methods_round2_d2 <- c(
  "WH-Empirical",
  "D1-Gaus-HC",
  "D1-Gaus-LC",
  "CD-0",
  "Naive"
)
table_round2_d2 <- out_table %>%
  dplyr::filter(`True Model` == "D2") %>%
  dplyr::select(
    `Ordinal Composite Type`,
    `True Settings`,
    `True Correlation`,
    `True Model`,
    `Component Prev Type`,
    dplyr::all_of(methods_round2_d2)
  )

for (comp_type in unname(comp_labels)) {
  
  file_stub <- paste0(
    "round2_d2_",
    stringr::str_to_lower(stringr::str_replace_all(comp_type, " ", "_"))
  )
  
  readr::write_csv(
    table_round2_d2 %>%
      dplyr::filter(`Ordinal Composite Type` == comp_type),
    file.path(table_dir, paste0(file_stub, ".csv"))
  )
  
  make_power_plot(
    data = out_long,
    comp_type = comp_type,
    true_model = "D2",
    methods_to_keep = methods_round2_d2,
    file_name = paste0(file_stub, ".png")
  )
}


# -----------------------------------------------------------------------------
# Print summary
# -----------------------------------------------------------------------------

message("Finished generating tables and figures.")
message("Tables saved to: ", table_dir)
message("Figures saved to: ", figure_dir)
message("Input file: ", summary_output_file)
message("Target n: ", target_n)

