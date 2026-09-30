#!/bin/bash
#SBATCH --job-name=hbcd-hdp
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/hbcd-hdp-%j.out
#SBATCH --error=logs/hbcd-hdp-%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs
module load python/3.11
module load R/4.5.1
make all
