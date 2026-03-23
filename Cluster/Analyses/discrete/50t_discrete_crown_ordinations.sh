#!/bin/bash
#SBATCH --job-name=50t_discrete_ordinations
#SBATCH --output=/users/bip24cns/acedisparity/discrete_crown/logs/50t_ord_%A_%a.out
#SBATCH --error=/users/bip24cns/acedisparity/discrete_crown/logs/50t_ord_%A_%a.err
#SBATCH --array=1-100          # 100 replicates
#SBATCH --nodes=1                     # 1 node per task
#SBATCH --ntasks=1                    # 1 R process per task
#SBATCH --cpus-per-task=1             # adjust if using parallel inside R
#SBATCH --mem=3G                     # adjust based on memory needs
#SBATCH --time=1:00:00               # 2 hours, adjust as needed
#SBATCH --mail-user=cnscutt1@sheffield.ac.uk
#SBATCH --mail-type=ALL

module load R/4.4.1-foss-2022b

export R_LIBS_USER=/users/$USER/R/x86_64-pc-linux-gnu-library/4.4

Rscript /users/$USER/acedisparity/discrete_crown/scripts/dynamic_discrete_crown_ordinations.R $SLURM_ARRAY_TASK_ID "50t" "8558401"