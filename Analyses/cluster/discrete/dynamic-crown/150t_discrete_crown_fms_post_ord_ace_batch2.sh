#!/bin/bash
#SBATCH --job-name=150t_post_ord_ace_batch_2
#SBATCH --output=/users/bip24cns/acedisparity/discrete_crown/logs/150t_post_ord_ace_batch2_%A_%a.out
#SBATCH --error=/users/bip24cns/acedisparity/discrete_crown/logs/150t_post_ord_ace_batch2_%A_%a.err
#SBATCH --array=1-500%50           # First 1000 tasks
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=96:00:00              # 96 hours for "all" levels
#SBATCH --mail-user=cnscutt1@sheffield.ac.uk
#SBATCH --mail-type=END,FAIL

module load R/4.4.1-foss-2022b
export R_LIBS_USER=/users/$USER/R/x86_64-pc-linux-gnu-library/4.4

Rscript /users/$USER/acedisparity/discrete_crown/scripts/dynamic_crown_post_ord_ace.R ${SLURM_ARRAY_TASK_ID} "150t" "8561712" "1000"