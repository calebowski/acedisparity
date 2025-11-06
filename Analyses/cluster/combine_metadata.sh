#!/bin/bash
#SBATCH --job-name=combine_metadata
#SBATCH --output=/users/$USER/acedisparity/trees/logs/combine_metadata_%A.out
#SBATCH --error=/users/$USER/acedisparity/trees/logs/combine_metadata_%A.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00
#SBATCH --mem=1G
#SBATCH --mail-user=cnscutt1@sheffield.ac.uk
#SBATCH --mail-type=ALL

module load R/4.4.1-foss-2022b
export R_LIBS_USER=/users/$USER/R/x86_64-pc-linux-gnu-library/4.4

Rscript /users/bip24cns/acedisparity/trees/scripts/combine_csv_metadata.R 1