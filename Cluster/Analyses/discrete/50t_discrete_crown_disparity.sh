#!/bin/bash
#SBATCH --job-name=50t_discrete_disparity
#SBATCH --output=/users/bip24cns/acedisparity/discrete_crown/logs/50t_discrete_disparity_%A_%a.out
#SBATCH --error=/users/bip24cns/acedisparity/discrete_crown/logs/50t_discrete_disparity_%A_%a.err
#SBATCH --array=1-100          
#SBATCH --nodes=1                     
#SBATCH --ntasks=1                   
#SBATCH --cpus-per-task=1              
#SBATCH --mem=8G                     
#SBATCH --time=3:00:00               
#SBATCH --mail-user=cnscutt1@sheffield.ac.uk
#SBATCH --mail-type=ALL

module load R/4.4.1-foss-2022b

export R_LIBS_USER=/users/$USER/R/x86_64-pc-linux-gnu-library/4.4


Rscript /users/$USER/acedisparity/discrete_crown/scripts/dynamic_crown_discrete_disparity.R $SLURM_ARRAY_TASK_ID "50t" "8558401"