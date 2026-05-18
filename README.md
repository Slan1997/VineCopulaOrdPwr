# VineCopulaOrdPwr

This repository contains R code for the simulation study and power analysis workflow accompanying the manuscript:

**Power and Sample Size Considerations for Ordinal Composite Endpoints in Clinical Trials**

The workflow supports:

- simulating correlated binary component outcomes under C-vine and D-vine copula structures;
- constructing ordinal composite outcomes from simulated binary components;
- estimating empirical ordinal composite PMFs and treatment-effect parameters;
- deriving vine-copula-based ordinal PMFs under different assumed copula structures;
- calculating Whitehead power using empirical, vine-copula-based, independence, and naive ordinal PMFs;
- generating manuscript-ready power tables and figures.

## Repository structure

```text
VineCopulaOrdPwr/
├── R/
│   ├── 00_settings_simulation.R
│   ├── 01_functions_simulation.R
│   ├── 02_run_simulation_array.R
│   ├── 02_run_simulation_array.slurm
│   ├── 03_functions_vco_pmf_power.R
│   ├── 04_combine_power_results.R
│   └── 05_generate_tables_figures.R
├── simulated_data/
├── figures/
├── tables/
├── .gitignore
└── VineCopulaOrdPwr.Rproj
