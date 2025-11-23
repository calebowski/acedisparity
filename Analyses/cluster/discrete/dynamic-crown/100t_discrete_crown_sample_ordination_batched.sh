#!/bin/bash
#SBATCH --job-name=100t_discrete_ordinations_batched
#SBATCH --output=/users/bip24cns/acedisparity/discrete_crown/logs/100t_discrete_sample_ord_batched_%A_%a.out
#SBATCH --error=/users/bip24cns/acedisparity/discrete_crown/logs/100t_discrete_sample_ord_batched_%A_%a.err
#SBATCH --array=1-1000%50  # 1,000 tasks, max 100 at once
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G            # 
#SBATCH --time=4:00:00      # 
#SBATCH --mail-user=cnscutt1@sheffield.ac.uk
#SBATCH --mail-type=ALL

module load R/4.4.1-foss-2022b
export R_LIBS_USER=/users/$USER/R/x86_64-pc-linux-gnu-library/4.4

Rscript /users/$USER/acedisparity/discrete_crown/scripts/dynamic_crown_sample_ordination_rate_batched.R ${SLURM_ARRAY_TASK_ID} "100t" "8561556"