#!/bin/bash
#SBATCH --job-name=lmm_combined_raw_rmaxes
#SBATCH --output=/users/bip24cns/acedisparity/discrete_crown/logs/lmm_combined_raw_rmaxes.out
#SBATCH --error=/users/bip24cns/acedisparity/discrete_crown/logs/lmm_combined_raw_rmaxes.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00
#SBATCH --mem=210G
#SBATCH --mail-user=cnscutt1@sheffield.ac.uk
#SBATCH --mail-type=ALL

# 1. Load the module
module load R/4.4.1-foss-2022b

# 2. Set the library path
export R_LIBS_USER=/users/bip24cns/R/x86_64-pc-linux-gnu-library/4.4

# 3. Run the test script
Rscript /users/bip24cns/acedisparity/discrete_crown/scripts/discrete_lmm.R "combined_raw_rmaxes"