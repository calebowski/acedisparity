#!/bin/bash
#SBATCH --job-name=acecontinuous100tord
#SBATCH --output=/users/bip24cns/acedisparity/continuous/logs/ord_continuous_100t_%A_%a.out
#SBATCH --error=/users/bip24cns/acedisparity/continuous/logs/ord_continuous_100t_%A_%a.err
#SBATCH --array=1-100       
#SBATCH --nodes=1                     # 1 node per task
#SBATCH --ntasks=1                    # 1 R process per task
#SBATCH --cpus-per-task=25              # adjust if using parallel inside R
#SBATCH --mem=24G                    # adjust based on memory needs
#SBATCH --time=80:00:00               
#SBATCH --mail-user=cnscutt1@sheffield.ac.uk
#SBATCH --mail-type=ALL

module load R/4.4.1-foss-2022b

export R_LIBS_USER=/users/$USER/R/x86_64-pc-linux-gnu-library/4.4


REPLICATE=${SLURM_ARRAY_TASK_ID}
MODELS=("bm" "bm_t" "ou_st" "ou_w" "ou_sh")
CORES_PER_MODEL=5
JOB_ID=8469104

echo "==== Starting replicate ${REPLICATE} on $(hostname) with ${SLURM_CPUS_PER_TASK} cores ===="
date


# Rscript /users/$USER/acedisparity/continuous/scripts/ace_continuous_100t_generate_tree.R ${REPLICATE}

# Run each model in parallel (background) and limit CPU usage per model
for MODEL in "${MODELS[@]}"; do
    echo "Starting model ${MODEL} for replicate ${REPLICATE} at $(date)"

    (
      export OMP_NUM_THREADS=$CORES_PER_MODEL
      Rscript /users/$USER/acedisparity/continuous/scripts/ace_continuous_process_model.R \
          "${REPLICATE}" "${MODEL}" "100t" ${JOB_ID} 
    ) &
    sleep 2
done

wait
echo "All models finished for replicate ${REPLICATE} at $(date)"
echo "Combining results..."
Rscript /users/$USER/acedisparity/continuous/scripts/combine_save_continuous_models.R ${REPLICATE} "100t" ${JOB_ID}

echo "Completed replicate ${REPLICATE}"

