#!/bin/bash
#SBATCH --job-name=test_for_loop_100t
#SBATCH --output=/users/bip24cns/acedisparity/continuous/logs/test_for_loop_%A_%a.out
#SBATCH --error=/users/bip24cns/acedisparity/continuous/logs/test_for_loop_%A_%a.err
#SBATCH --array=1-100       # Test just first 5 replicates
#SBATCH --nodes=1                     
#SBATCH --ntasks=1                    
#SBATCH --cpus-per-task=25              
#SBATCH --mem=24G                    
#SBATCH --time=00:10:00               
#SBATCH --mail-user=cnscutt1@sheffield.ac.uk
#SBATCH --mail-type=ALL

module load R/4.4.1-foss-2022b
export R_LIBS_USER=/users/$USER/R/x86_64-pc-linux-gnu-library/4.4

REPLICATE=${SLURM_ARRAY_TASK_ID}
MODELS=("bm" "bm_t" "ou_st" "ou_w" "ou_sh")
CORES_PER_MODEL=5

echo "==== Starting replicate ${REPLICATE} on $(hostname) with ${SLURM_CPUS_PER_TASK} cores ===="
date

# Replicate your exact for loop structure
for MODEL in "${MODELS[@]}"; do
    echo "Starting model ${MODEL} for replicate ${REPLICATE} at $(date)"

    # Run in background with the same structure as your main script
    (
      export OMP_NUM_THREADS=$CORES_PER_MODEL
      # Use the debug script instead of the full process script
      Rscript /users/$USER/acedisparity/continuous/scripts/read_tree_debug.R \
          "${REPLICATE}" "${MODEL}" "100t"
    ) &

    sleep 2
done

wait
echo "All models finished for replicate ${REPLICATE} at $(date)"
echo "Completed replicate ${REPLICATE}"                                                                                                                                                                                                                                                                                                                                                                                                                                                                    