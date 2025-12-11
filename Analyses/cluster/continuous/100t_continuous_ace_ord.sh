#!/bin/bash
#SBATCH --job-name=100t_continuous_ace_ord
#SBATCH --output=/users/bip24cns/acedisparity/continuous/logs/100t_continuous_ace_ord_%A_%a.out
#SBATCH --error=/users/bip24cns/acedisparity/continuous/logs/100t_continuous_ace_ord_%A_%a.err
#SBATCH --array=1-100      
#SBATCH --nodes=1                     # 1 node per task
#SBATCH --ntasks=1                    # 1 R process per task
#SBATCH --cpus-per-task=1            # adjust if using parallel inside R
#SBATCH --mem=12G                    # adjust based on memory needs
#SBATCH --time=24:00:00               
#SBATCH --mail-user=cnscutt1@sheffield.ac.uk
#SBATCH --mail-type=ALL

module load R/4.4.1-foss-2022b

export R_LIBS_USER=/users/$USER/R/x86_64-pc-linux-gnu-library/4.4

Rscript /users/$USER/acedisparity/continuous/scripts/dynamic_continuous_crown_pre_ord.R $SLURM_ARRAY_TASK_ID "100t"