#!/bin/bash
#SBATCH --job-name=tree_creation_100
#SBATCH --output=/users/bip24cns/acedisparity/continuous/logs/tree_continuous_100t_%A_%a.out
#SBATCH --error=/users/bip24cns/acedisparity/continuous/logs/tree_continuous_100t_%A_%a.err
#SBATCH --array=1-100       
#SBATCH --nodes=1                     # 1 node per task
#SBATCH --ntasks=1                    # 1 R process per task
#SBATCH --cpus-per-task=1              # adjust if using parallel inside R
#SBATCH --mem=1G                    # adjust based on memory needs
#SBATCH --time=00:20:00               
#SBATCH --mail-user=cnscutt1@sheffield.ac.uk
#SBATCH --mail-type=ALL

module load R/4.4.1-foss-2022b

export R_LIBS_USER=/users/$USER/R/x86_64-pc-linux-gnu-library/4.4
REPLICATE=${SLURM_ARRAY_TASK_ID}
Rscript /users/$USER/acedisparity/continuous/scripts/ace_continuous_100t_generate_tree.R ${REPLICATE}