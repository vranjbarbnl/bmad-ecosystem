#!/usr/bin/env bash
#SBATCH --job-name=fccee_z_dm_chi002
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

# Add SDCC-approved --account, --partition, --time, and memory directives before submission.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit this script with sbatch from the run directory}"
./run.sh
