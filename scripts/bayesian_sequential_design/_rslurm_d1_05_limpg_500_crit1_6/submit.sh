#!/bin/bash
#
#SBATCH --array=0-0
#SBATCH --cpus-per-task=16
#SBATCH --job-name=d1_05_limpg_500_crit1_6
#SBATCH --output=slurm_%a.out
/usr/lib64/R/bin/Rscript --vanilla slurm_run.R
