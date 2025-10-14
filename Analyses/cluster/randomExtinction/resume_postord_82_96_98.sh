#!/bin/bash
#SBATCH --job-name=resume_postord_82_96_98
#SBATCH --output=/users/bip24cns/acedisparity/randomExtinction/logs/resume_postord_82_96_98_%A_%a.out
#SBATCH --error=/users/bip24cns/acedisparity/randomExtinction/logs/resume_postord_82_96_98_%A_%a.err
#SBATCH --array=82,96,98
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=64G
#SBATCH --time=80:00:00
#SBATCH --mail-user=cnscutt1@sheffield.ac.uk
#SBATCH --mail-type=ALL

module load R/4.4.1-foss-2022b
export R_LIBS_USER=/users/$USER/R/x86_64-pc-linux-gnu-library/4.4

Rscript /users/$USER/acedisparity/randomExtinction/scripts/resume_150t_replicates.R $SLURM_ARRAY_TASK_ID