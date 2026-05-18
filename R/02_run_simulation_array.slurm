#!/bin/bash
#SBATCH --job-name=sim_array_1e5
#SBATCH --time=5-22:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --array=1-30
#SBATCH --output=slurm_out/slurm-%A_%a.out

# -----------------------------------------------------------------------------
# 02_run_simulation_array.slurm
# -----------------------------------------------------------------------------
# Run the large-n simulation array.
#
# Each array task corresponds to one row of the scenario grid defined in
# 00_settings_simulation.R. The task index is passed through SLURM_ARRAY_TASK_ID
# and used by 02_run_simulation_array.R.
#
# Optional email notifications:
#   #SBATCH --mail-user=your.email@example.com
#   #SBATCH --mail-type=ALL
# -----------------------------------------------------------------------------

mkdir -p slurm_out

setup_accre_software_stack
module load gcc r

Rscript --no-save --no-restore 02_run_simulation_array.R