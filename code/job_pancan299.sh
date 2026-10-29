#!/bin/bash
#SBATCH --job-name=pancan299
#SBATCH --output=/dartfs-hpc/rc/home/m/f0052zm/szhao_lab_share/diffDriver/data_run_pancancer299/logs/t_%A_%a.out
#SBATCH --error=/dartfs-hpc/rc/home/m/f0052zm/szhao_lab_share/diffDriver/data_run_pancancer299/logs/t_%A_%a.err
#SBATCH --array=1-237%25
#SBATCH --time=24:00:00
#SBATCH --mem=72G
#SBATCH --cpus-per-task=4
set -euo pipefail
module load R
cd /dartfs-hpc/rc/home/m/f0052zm/szhao_lab_share/diffDriver/data_run_pancancer299
Rscript run_cohort_multi.R
