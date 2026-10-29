#!/bin/bash
#SBATCH --job-name=pancan299_rec
#SBATCH --output=/dartfs-hpc/rc/home/m/f0052zm/szhao_lab_share/diffDriver/data_run_pancancer299/logs/r_%A_%a.out
#SBATCH --error=/dartfs-hpc/rc/home/m/f0052zm/szhao_lab_share/diffDriver/data_run_pancancer299/logs/r_%A_%a.err
#SBATCH --array=1-32%16
#SBATCH --time=48:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=4
set -euo pipefail
module load R
cd /dartfs-hpc/rc/home/m/f0052zm/szhao_lab_share/diffDriver/data_run_pancancer299
Rscript recover_cohort.R
