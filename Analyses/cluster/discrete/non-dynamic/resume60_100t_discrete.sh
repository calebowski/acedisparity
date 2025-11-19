#!/bin/bash
#SBATCH --job-name=resume_rep_60_100t
#SBATCH --output=/users/bip24cns/acedisparity/discrete/logs/resume_rep_60_100t.out
#SBATCH --error=/users/bip24cns/acedisparity/discrete/logs/resume_rep_60_100t.err
#SBATCH --nodes=1                     # 1 node per task
#SBATCH --ntasks=1                    # 1 R process per task
#SBATCH --cpus-per-task=30             # adjust if using parallel inside R
#SBATCH --mem=12G                     # adjust based on memory needs
#SBATCH --time=85:00:00               # 8 hours, adjust as needed
#SBATCH --mail-user=cnscutt1@sheffield.ac.uk
#SBATCH --mail-type=ALL

module load R/4.4.1-foss-2022b

export R_LIBS_USER=/users/$USER/R/x86_64-pc-linux-gnu-library/4.4

Rscript /users/$USER/acedisparity/discrete/scripts/60_finish_job.R 1