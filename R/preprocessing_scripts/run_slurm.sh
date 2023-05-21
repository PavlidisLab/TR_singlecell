#!/bin/bash
#SBATCH --job-name=run_r_scripts
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=365-0:0:0   # 365 days (1 year) time limit

# Load R module
module load R

# Set paths
SCRIPTS_DIR=/home/amorin/Projects/TR_singlecell/R/test_slurm/

# Loop through R scripts and submit jobs
for SCRIPT in "$SCRIPTS_DIR"/*.R; do
  Rscript --vanilla $SCRIPT
done
